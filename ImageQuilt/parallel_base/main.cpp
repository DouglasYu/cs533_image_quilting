#include <iostream>
#include "ImageQuilt.h"

int main()
{
	ImageQuilt* image_quilt = new ImageQuilt("../input/text.bmp", "../output/", 64, 5, 16, 0.1);
	image_quilt->synthesize();
	delete image_quilt;
	return 0;
}