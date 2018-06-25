
#include "r3d.h"
#include <stdio.h>
#include <time.h>

//The amount of times the program will be performed for timer.
#define REP_TIMES 1000000

// A simple printer that prints all BREP information provided to it
void printer(int faces, int* verts, int** faceindices)
{
	printf("Number of faces: %d\n", faces);
	printf("# of vertices by face: \n");	
	for (int k = 0; k<faces; ++k)
		printf("%d ", *(verts+k));
	printf("\nVertex indices: \n");	
	for (int i = 0; i<faces; ++i)
	{
		for (int j = 0; j<*(verts+i); ++j)
		{
			printf("%d ", *(*(faceindices+i)+j));
		}
		printf("\n");
	}
}

//This is for a cube!
int main(){

	int nfaces = 6; // number of faces in poly
	int nvertsperface[] = {4,4,4,4,4,4}; // number of vertices per face


	r3d_rvec3 verts[] = {
		{0,0,0},{1,0,0},{1,1,0},{0,1,0},
		{0,0,1},{1,0,1},{1,1,1},{0,1,1}}; //the vertices
	int nverts = sizeof(verts)/sizeof(verts)[0];
	
	// the faces, holds the brep of each face
	int f0[] = {2,6,5,1};
	int f1[] = {1,5,4,0};
	int f2[] = {4,7,3,0};
	int f3[] = {7,6,2,3};
	int f4[] = {6,7,4,5};
	int f5[] = {3,2,1,0};
	int* faceinds[] = {f0,f1,f2,f3,f4,f5}; // the list of all faces as defined by the brep

	r3d_poly cube; // the cube
	r3d_plane planes[] = {{{1,0,0}, -0.5}}; // the planes
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes

	// variables for clipping and converting to brep
	r3d_brep *cuberep;
	r3d_int numcomponents;
	
	//timer
    clock_t start = clock();
	for (int i = 0; i<REP_TIMES; i++)
		{
		r3d_init_poly(&cube, verts, nverts, faceinds, nvertsperface, nfaces); 
		r3d_clip(&cube, planes, nplanes);
		r3d_init_brep(&cube, &cuberep, &numcomponents);
		}
	clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed);
	
	/*
	r3d_init_poly(&cube, verts, nverts, faceinds, nvertsperface, nfaces);
	r3d_print(&cube);
	printf("\n");
	r3d_clip(&cube, planes, nplanes);
	r3d_print(&cube);
	printf("\n");
	r3d_print_brep(&cuberep, numcomponents);
	*/
	
}
