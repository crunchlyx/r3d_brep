#include "r2d.h"
#include <stdio.h>
#include <time.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 10000000

//This is for a square!
int main(){

	r2d_rvec2 verts[] = {{0,0},{1,0},{1,1},{0,1}}; //the vertices in counterclockwise order
	int nverts = sizeof(verts)/sizeof(verts)[0];

	r2d_poly square; // the square
	r2d_plane planes[] = {{{1,0}, -0.5}}; // the planes
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes
	
	//timer
    clock_t start = clock();
	for (int i = 0; i<REP_TIMES; i++)
		{
		r2d_init_poly(&square, verts, nverts); 
		r2d_clip(&square, planes, nplanes);
		}
	clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed);	
	/*
	r2d_init_poly(&square, verts, nverts); 
	r2d_print(&square);
	printf("\n");
	r2d_clip(&square, planes, nplanes);
	r2d_print(&square);
	*/
}
