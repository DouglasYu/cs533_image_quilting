#include "ImageQuilt.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include <queue>
#include <stdint.h>
#define STBI_ONLY_BMP
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define USE_CIELAB

#include <functional>

using std::cout;
using std::endl;
using std::vector;
using std::pair;

uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

template <class Tp>
inline void DoNotOptimize(Tp const& value) {
  asm volatile("" : : "r,m"(value) : "memory");
}

template <class Tp>
inline void DoNotOptimize(Tp& value) {
#if defined(__clang__)
  asm volatile("" : "+r,m"(value) : : "memory");
#else
  asm volatile("" : "+m,r"(value) : : "memory");
#endif
}

uint64_t global_start_time;


ImageQuilt::~ImageQuilt()
{
}

inline void ImageQuilt::fillPatch(Patch * patch, unsigned int inner_w, unsigned int inner_h, unsigned int out_w, unsigned int out_h)
{
	patch->image->l(inner_w, inner_h, output_image->l(out_w, out_h));
	patch->image->a(inner_w, inner_h, output_image->a(out_w, out_h));
	patch->image->b(inner_w, inner_h, output_image->b(out_w, out_h));
}

inline void ImageQuilt::applyPatch(Patch * patch, unsigned int out_w, unsigned int out_h, unsigned int inner_w, unsigned int inner_h)
{
	output_image->l(out_w, out_h, patch->image->l(inner_w, inner_h));
	output_image->a(out_w, out_h, patch->image->a(inner_w, inner_h));
	output_image->b(out_w, out_h, patch->image->b(inner_w, inner_h));
}

// Added
void ImageQuilt::put_first_tile() {
	printf("Tile:  0, 0 |  1/64 | thread  x\n");
	isDone[0][0] = rdtsc() - global_start_time;

	// pick a random patch, apply selected patch to output starting at [0,0]
	Patch* initial = randomPatch();
	input_image->rgb2xyz();
	input_image->xyz2cielab();
	output_image->rgb2xyz();
	output_image->xyz2cielab();

	unsigned int line = 0;
	unsigned int col = 0;
	for (unsigned int i = initial->hLocation; i < initial->hLocation + tilesize; i++) {
		for (unsigned int j = initial->wLocation; j < initial->wLocation + tilesize; j++) {
			output_image->l(col, line, input_image->l(j, i));
			output_image->a(col, line, input_image->a(j, i));
			output_image->b(col, line, input_image->b(j, i));
			col++;
		}
		line++;
		col = 0;
	}
}

// Added
void ImageQuilt::put_tile(unsigned int tile_hi, unsigned int tile_wi, int id) {
	// coordinates on the output image
	auto output_w = tile_wi*tilesize - tile_wi*overlap;
	auto output_h = tile_hi*tilesize - tile_hi*overlap;

	/** 
	 * 1. scan input image 
	 * For every location, search the input texture for a set of blocks that satisfy the overlap constraints(above and left) within some error tolerance.
	 * Randomly pick one such block.
	 */
	double lowest = DBL_MAX;
	vector<pair<pair<int, int>, double>> distances;
	uint64_t inner_loop_count = 0;
	// bool found_zero = false;
	for (unsigned int img_h = 0; img_h < input_image->height - tilesize; img_h++) {
	for (unsigned int img_w = 0; img_w < input_image->width - tilesize; img_w++) {
		inner_loop_count++;
		/**
		 * 1A. Find a region in the input picture that have the lowest error ratio with the 
		 * current output picture overlap portion
		 */
		double sum = 0.0;
		if(tile_hi == 0 && tile_wi > 0) //First Row
		{
			for (unsigned int inner_h = 0; inner_h < tilesize; inner_h++) {
			for (unsigned int inner_w = 0; inner_w < overlap; inner_w++) {
				sum += errorCalc(img_w + inner_w, img_h + inner_h, output_w + inner_w, output_h + inner_h);
			}}
		}
		else if(tile_hi > 0 && tile_wi == 0) //First Column
		{
			for (unsigned int inner_h = 0; inner_h < overlap; inner_h++) {
			for (unsigned int inner_w = 0; inner_w < tilesize; inner_w++) {
				sum += errorCalc(img_w + inner_w, img_h + inner_h, output_w + inner_w, output_h + inner_h);
			}}
		}
		else //Normal situations
		{
			for (unsigned int inner_h = 0; inner_h < overlap; inner_h++) {
			for (unsigned int inner_w = 0; inner_w < tilesize; inner_w++) {
				sum += errorCalc(img_w + inner_w, img_h + inner_h, output_w + inner_w, output_h + inner_h);
			}}
			for (unsigned int inner_h = overlap; inner_h < tilesize; inner_h++) {
			for (unsigned int inner_w = 0; inner_w < overlap; inner_w++) {
				sum += errorCalc(img_w + inner_w, img_h + inner_h, output_w + inner_w, output_h + inner_h);
			}}
		}

		sum = sqrt(sum);
		if(sum == 0.0)	//Perfect match
		{
			distances.clear();
			pair<int, int> location(img_w, img_h);
			pair<pair<int, int>, double> element(location, sum);
			distances.push_back(element);
			// found_zero = true;
			// cout << "zero" << endl;
			goto scan_finish;
		}
		else if (sum <= (1+error)*lowest) {
			if (sum < lowest) {
				lowest = sum;
				// cout << "Lowest: " << lowest << endl;
			}
			pair<int, int> location(img_w, img_h);
			pair<pair<int, int>, double> element(location, sum);
			distances.push_back(element);
		}
	}}

	scan_finish:
	// count how many patches are below the cutoff
	auto cutoff = (1 + error)*lowest;
	auto matches = -1;
	for (vector<pair<pair<int, int>, double>>::reverse_iterator it = distances.rbegin(); it != distances.rend(); ++it) {
		if (it->second > cutoff) break;
		matches++;
		inner_loop_count++;
	}

	/** 
	 * 2. From all the potential choice candidates[Ones below the err cutoff]
	 * we randomly choose one to be the resulting tile
	 */
	unsigned int randNum = rand() % (matches - 0 + 1);
	//cout << "Picked: " << randNum << endl;
	unsigned int patch_w = distances[distances.size() - randNum - 1].first.first;
	unsigned int patch_h = distances[distances.size() - randNum - 1].first.second;
	double err = distances[distances.size() - randNum - 1].second;
	// printf("%d/%d (thread %d) ", tile_hi*num_tiles + tile_wi + 1, total_tiles, id);
	// cout << id << "--Picked: " << randNum << " out of: " << matches << " with error: " << err << " lowest: " << lowest << " inner_loop_count: " << inner_loop_count << endl;
	// printf("[Test] Thread %d completed (%d, %d)\n", id, tile_wi, tile_hi);
	
	/** 
	 * 3. Create a tempoary tile[Patch] as the workstation.
	 * On this workstation, we find out a cutline in the overlay portion, outside this cutline, we place the existing portion 
	 * of the output image; within this cutline, we place the portion from the input_image.
	 */
	Patch* patch = new Patch(patch_w, patch_h, tilesize, tilesize, input_image);
	patch->image->rgb2xyz();
	patch->image->xyz2cielab();

	/* Only do patching if we have err, else we could directly copy and paste */
	if (err > 0) {

		if (tile_wi > 0) //left overlap
		{
			auto final_cut = minCut(output_w, output_h, patch_w, patch_h, true);
			for (unsigned int inner_h = 0; inner_h < tilesize; inner_h++){
			for (int inner_w = 0; inner_w < final_cut->at(inner_h); inner_w++){
				fillPatch(patch, inner_w, inner_h, output_w + inner_w, output_h + inner_h);
			}}
			delete final_cut;
		}
		
		if (tile_hi > 0) // top overlap
		{
			auto final_cut = minCut(output_w, output_h, patch_w, patch_h, false);
			for (unsigned int inner_w = 0; inner_w < tilesize; inner_w++){
			for (int inner_h = 0; inner_h < final_cut->at(inner_w); inner_h++){
				fillPatch(patch, inner_w, inner_h, output_w + inner_w, output_h + inner_h);
			}}
			delete final_cut;
		}
	}

	/** 
	 * 4. We copy the content from this tempoary workstation[Patch] to the output image
	 */
	for (unsigned int inner_h = 0; inner_h < tilesize; inner_h++){
	for (unsigned int inner_w = 0; inner_w < tilesize; inner_w++){
		applyPatch(patch, output_w + inner_w, output_h + inner_h, inner_w, inner_h);
	}}

	delete patch;
}

// Added
bool ImageQuilt::updateCursor(int *x, int *y, int offset) {
	int xn = *x - offset;
	int yn = *y + offset;
	int itcp = xn + yn;
	int delta;
	while (xn < 0 || yn >= num_tiles) {
		if (itcp < num_tiles - 1) {
			itcp += 1;
			delta = - xn - 1;
			xn = itcp - delta;
			yn = delta;
		} else if (itcp >= num_tiles - 1 && itcp < 2 * num_tiles - 2) {
			itcp += 1;
			delta = yn - num_tiles;
			xn = num_tiles - 1 - delta;
			yn = itcp - num_tiles + delta + 1;
		} else {
			return false;
		}
	}
	// printf("[Test] (%d, %d)/%d + %d = (%d, %d)/%d\n", *x, *y, *y*num_tiles+*x, offset, xn, yn, yn*num_tiles+xn);
	*x = xn;
	*y = yn;
	return true;
}

// void ImageQuilt::printMtx() {
// 	for (int i = 0; i < num_tiles; i++) {
// 		for (int j = 0; j < num_tiles; j++) {
// 			if (isDone[i][j] == true)
// 				printf("1 ");
// 			else
// 				printf("0 ");
// 		}
// 		printf("\n");
// 	}
// }


/**
 * This is the pre-define location method..
*/
void ImageQuilt::put_tile_thread (int id) {
	int tile_wi = 0;	// x
	int tile_hi = 0; 	// y
	uint64_t stall_count = 0;	//Not used, BUT CAN'T BE DELETED!!!
	uint64_t total_time = 0;
	uint64_t inefficient_time = 0;

	// std::vector<int> v;

	updateCursor(&tile_wi, &tile_hi, id+1);
	do {

		uint64_t start_time = rdtsc();
		while((tile_wi != 0) && (0 == (isDone[tile_hi][tile_wi-1]))) DoNotOptimize(stall_count++);	//DO NOT DELETE, serve as volatile
		while((tile_hi != 0) && (0 == (isDone[tile_hi-1][tile_wi]))) DoNotOptimize(stall_count++);
		uint64_t time_ms = (rdtsc() - start_time)/1000000/2.6;


		// printMtx();
		uint64_t exe_begin_time = rdtsc();
		put_tile((unsigned int)tile_hi, (unsigned int)tile_wi, id);
		uint64_t exe_end_time = (rdtsc() - exe_begin_time)/1000000/2.6;
		printf("Tile: %2d,%2d | %2d/%2d | thread %2d | execution_time_ms %llu | stall_time_ms %llu | curr_time_ms %llu\n", 
			tile_wi, tile_hi, tile_hi*num_tiles + tile_wi + 1, total_tiles, id, exe_end_time,time_ms, (uint64_t)((rdtsc() - global_start_time)/1000000/2.6));
		isDone[tile_hi][tile_wi] = rdtsc() - global_start_time;
		stall_count = 0;

		total_time += (time_ms + exe_end_time);

	} while (updateCursor(&tile_wi, &tile_hi, num_threads));

	mtx.lock();
	std::cout << "ID -- " << id << " Total Time -- " << total_time/1000 << std::endl;
	mtx.unlock();
}


/**
 * This is the dynamic queue method
 * AKA we use lock free atomic operation queue to 
*/



// Modified
void ImageQuilt::synthesize() {
	put_first_tile();
	global_start_time = rdtsc();

	std::vector<std::thread> threads;
	for (int i = 0; i < num_threads; i++) {
        threads.emplace_back([this, i] {put_tile_thread(i);});
    }

	for (auto& thread : threads) {
        thread.join();
    }

	// store the result to .bmp
	output_image->cielab2xyz();
	output_image->xyz2rgb();
	stbi_write_bmp(output_filename.c_str(), output_width, output_height, 3, output_image->data);

	// cleanup
	delete input_image;
	delete output_image;
}

Image* ImageQuilt::get_output() const
{
	return output_image;
}

void ImageQuilt::loadImage()
{
	int x, y, n;
	uint8_t* data = stbi_load(input_filename.c_str(), &x, &y, &n, 0);
	assert(data);
	input_image = new Image(x, y, data);
}

Patch * ImageQuilt::randomPatch()
{
	unsigned int w = 0;//rand() % (input_image->width - tilesize);
	unsigned int h = 0;//rand() % (input_image->height - tilesize);
	Patch* p = new Patch(w, h, tilesize, tilesize, input_image);
	return p;
}

// using L2 Norm
inline double ImageQuilt::errorCalc(unsigned int in_w, unsigned int in_h, unsigned out_w, unsigned int out_h) const
{

	double val = 0.0;
	val += pow(input_image->l(in_w, in_h) - output_image->l(out_w, out_h), 2);
	val += pow(input_image->a(in_w, in_h) - output_image->a(out_w, out_h), 2);
	val += pow(input_image->b(in_w, in_h) - output_image->b(out_w, out_h), 2);
	return val;
}

vector<int>* ImageQuilt::minCut(unsigned int pos_w, unsigned int pos_h, unsigned int patch_w, unsigned int patch_h, const bool left)
{
	unsigned int vertical, horizontal;
	if(left)
	{
		vertical = tilesize;
		horizontal = overlap;
	}
	else
	{
		vertical = overlap;
		horizontal = tilesize;
	}

	int** dis = new int*[vertical];
	bool** visited = new bool*[vertical];
	double** overlap_diff = new double*[vertical];
	for (unsigned int h = 0; h < vertical; h++)
	{
		dis[h] = new int[horizontal];
		visited[h] = new bool[horizontal];
		overlap_diff[h] = new double[horizontal];
		for (unsigned int w = 0; w < horizontal; w++)
		{
			dis[h][w] = INT_MAX;
			visited[h][w] = false;
			/*double err = 0.3*pow(input_image->rc(patch_w + w, patch_h + h) - output_image->rc(pos_w + w, pos_h + h), 2);
			err += 0.59*pow(input_image->gc(patch_w + w, patch_h + h) - output_image->gc(pos_w + w, pos_h + h), 2);
			err += 0.11*pow(input_image->bc(patch_w + w, patch_h + h) - output_image->bc(pos_w + w, pos_h + h), 2);*/
			overlap_diff[h][w] = sqrt(errorCalc(patch_w + w, patch_h + h, pos_w + w, pos_h + h));
		}
	}
	typedef pair<pair<int, int>, vector<int>*> tt;
	std::priority_queue< pair<tt, double>, vector<pair<tt, double>>, Prioritize<tt> > pq;
	vector<int>* final_cut(0);
	if (left)
	{
		for (unsigned int w = 0; w < horizontal; w++)
		{
			pq.push(std::make_pair(std::make_pair(std::make_pair(w, 0), new vector<int>(1, w)), overlap_diff[0][w]));
		}
		while (!pq.empty())
		{
			auto current = pq.top(); //Current vertex. The shortest distance for this has been found
			pq.pop();
			unsigned int current_w = current.first.first.first;
			unsigned int current_h = current.first.first.second;
			vector<int>* current_vector = current.first.second;
			auto current_error = current.second;
			if (current_h == vertical - 1)
			{
				// found the best
				final_cut = current_vector;
				break;
			}
			if (visited[current_h][current_w])
			{
				delete current_vector;
				continue;
			}
			visited[current_h][current_w] = true;
			unsigned int start = (current_w - 1) >= 0 ? (current_w - 1) : 0;
			unsigned int end = (current_w + 1) < horizontal ? (current_w + 1) : current_w;
			for (unsigned int i = start;i <= end;i++)
			{
				if (!visited[current_h + 1][i] && overlap_diff[current_h + 1][i] + current_error < dis[current_h + 1][i])
				{
					auto entry = std::make_pair(std::make_pair(i, current_h + 1), new vector<int>(current_vector->begin(), current_vector->end()));
					entry.second->push_back(i);
					pq.push(std::make_pair(entry, dis[current_h + 1][i] = overlap_diff[current_h + 1][i] + current_error));
				}
			}
			delete current_vector;
		}
	}
	else
	{
		for (unsigned int h = 0; h < vertical; h++)
		{
			pq.push(std::make_pair(std::make_pair(std::make_pair(0, h), new vector<int>(1, h)), overlap_diff[h][0]));
		}
		while (!pq.empty())
		{
			auto current = pq.top(); //Current vertex. The shortest distance for this has been found
			pq.pop();
			unsigned int current_w = current.first.first.first;
			unsigned int current_h = current.first.first.second;
			vector<int>* current_vector = current.first.second;
			auto current_error = current.second;
			if (current_w == horizontal - 1)
			{
				// found the best
				final_cut = current_vector;
				break;
			}
			if (visited[current_h][current_w])
			{
				delete current_vector;
				continue;
			}
			visited[current_h][current_w] = true;
			unsigned int start = (current_h - 1) >= 0 ? (current_h - 1) : 0;
			unsigned int end = (current_h + 1) <= vertical ? (current_h + 1) : current_h;
			for (auto i = start;i < end;i++)
			{
				if (!visited[i][current_w + 1] && overlap_diff[i][current_w + 1] + current_error < dis[i][current_w + 1])
				{
					auto entry = std::make_pair(std::make_pair(current_w + 1, i), new vector<int>(current_vector->begin(), current_vector->end()));
					entry.second->push_back(i);
					pq.push(std::make_pair(entry, dis[i][current_w + 1] = overlap_diff[i][current_w + 1] + current_error));
				}
			}
			delete current_vector;
		}
	}
	// cleanup
	for (unsigned int h = 0; h < vertical; h++)
	{
		delete dis[h];
		delete visited[h];
		delete overlap_diff[h];
	}
	delete dis, visited, overlap_diff;
	return final_cut;
}