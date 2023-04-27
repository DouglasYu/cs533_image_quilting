#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include "custom_types.h"
#include <math.h>       /* round, floor, ceil, trunc */
#include <queue>
#include <utility>

// Added
#include <thread>
#include <mutex>
#include <barrier>

using std::string;
class ImageQuilt
{
public:
	ImageQuilt(const string& _input_filename, const string& _output_folder, unsigned _tilesize, unsigned _num_tiles, unsigned _overlap, double _error)
		: input_filename(_input_filename),
		  tilesize(_tilesize),
		  num_tiles(_num_tiles),
		  overlap(_overlap),
		  error(_error),
		  sync_point(std::thread::hardware_concurrency())

	{
		if (overlap == 0)
		{
			overlap = static_cast<unsigned int>(round(tilesize / 6));
		}
		if(error == 0)
		{
			error = 0.1;
		}
		output_filename = _output_folder + "output-" + std::to_string(tilesize) + "-" + std::to_string(num_tiles) + "-" + std::to_string(overlap) + ".bmp";
		output_width = num_tiles * tilesize - ((num_tiles - 1) * overlap);
		output_height = num_tiles * tilesize - ((num_tiles - 1) * overlap);
		output_image = new Image(output_width, output_height);

		//Add
		num_threads = std::thread::hardware_concurrency() - 1;	


		srand(time(NULL));
		loadImage();
	}
	~ImageQuilt();
	inline void fillPatch(Patch* patch, unsigned int inner_w, unsigned int inner_h, unsigned int out_w, unsigned int out_h);
	inline void applyPatch(Patch* patch, unsigned int out_w, unsigned int out_h, unsigned int inner_w, unsigned int inner_h);
	void synthesize();
	Image* get_output() const;
private:
	void loadImage();
	Patch* randomPatch();
	inline double errorCalc(unsigned int in_w, unsigned int in_h, unsigned out_w, unsigned int out_h) const;
	std::vector<int>* minCut(unsigned int pos_w, unsigned int pos_h, unsigned int patch_w, unsigned int patch_h, const bool left);
	void search_worker_thread (int id);

	//Thread
	int num_threads;		// Number of CPU cores available, also number of threads to launch

	//Synch related
	std::atomic<int32_t> perfect_match;	//used for optimization
	std::barrier<decltype([]() noexcept {})> sync_point;
	std::mutex mtx;
	std::vector<std::pair<std::pair<int, int>, double>> distances;


	string input_filename = {};
	string output_filename = {};
	unsigned int tilesize = {};
	unsigned int num_tiles = {};
	unsigned int output_width = {};
	unsigned int output_height = {};
	unsigned int overlap = {};
	double error = {};
	Image* input_image = {};
	Image* output_image = {};
};

