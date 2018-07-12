#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
//#define REP_TIMES 10000000
#define REP_TIMES 1
//#define REP_TIMES 1

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

// Buffer limit
#define MAX_VERTS 64
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

//free allocated memory
/*void free_poly(r2d_brep* poly) {
	r2d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		free(poly->verts);
	}
	free(poly->nverts);
	free(poly);
}*/

//clips a poly, takes a statically allocated buffer and works in place
void inplace_clip_poly(r2d_brep* poly, r2d_plane* planes, r2d_int nplanes) {

	const r2d_real ZERO = 0.0;	
	r2d_int vnext, onv, workindex, newindex;
	
	//loop over each clipping plane
	for(int p = 0; p < nplanes; p++) {

		r2d_rvec2 *verts = poly->verts;
		r2d_real sdists[onv];
		onv=workindex=poly->nverts;
		//calculate signed distances
		for(int v = 0; v < onv; ++v) {
			sdists[v] = (planes+p)->d + dot(verts[v], (planes+p)->n);
		}

		//loop over each vertex
		for(int v = 0; v < onv; ++v) {

			//note the indices of vertices that will be kept
			if (sdists[v] >= 0){
				verts[workindex++] = verts[v];
			}
			vnext = (v+1)%onv;
			
			//interpolate new vertex if edge crosses
			if (sdists[v] * sdists[vnext] < ZERO) {
				r2d_real alpha = sdists[vnext] / (sdists[vnext] - sdists[v]);
				verts[workindex].x = alpha*verts[v].x + (1-alpha)*verts[vnext].x;
				verts[workindex].y = alpha*verts[v].y + (1-alpha)*verts[vnext].y;
				workindex++;
			}	
		}

		newindex=0;
		for(int i = onv; i < workindex; ++i) {
			verts[newindex++] = verts[i];
		}
		poly->nverts = workindex-onv;
	}
}

//generate a random poly in a buffer of the right size, for use by clip_poly and r2d_clip
void init_random_poly(r2d_rvec2* verts, r2d_int* nverts) {

	//generate random number of vertices
	int randint = (rand()%8)+3;
	*nverts = randint;

	for(int i = 0; i < randint; ++i){
		r2d_real randx = RAND_COORD_LIMIT*(double)rand()/(double)RAND_MAX;
		r2d_real randy = RAND_COORD_LIMIT*(double)rand()/(double)RAND_MAX;
		r2d_rvec2 randvert;
		randvert.x = randx;
		randvert.y = randy;
		verts[i] = randvert;
	}
}

void info_poly2(r2d_brep* poly) { 
	printf("\n______________________________________________\n");
	
	r2d_rvec2* verts = poly->verts;
	r2d_int nverts = poly->nverts;
	//number of vertices	
	printf("%d vertices in poly: \n", nverts);
	
	//all vertex coordinates
	printf("Vertex coordinates:\n");
	for (int i = 0; i < nverts; ++i) {
		r2d_rvec2 printvert = verts[i];

		printf("%f ", printvert.x);
		printf("%f\n", printvert.y);
	}
	printf("______________________________________________\n");
}


int main(){

	//number of vertices
	r2d_int nverts;
	srand(time(0));
	r2d_rvec2 randverts[R2D_MAX_VERTS];
	// the intersection planes
	r2d_plane planes[] = {{{1,0}, -5.0}, {{0,1}, -5.0}};
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes
	// simple timer loop	
	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) 
		{
		
		// generate a random poly, used by all functions
		init_random_poly(randverts, &nverts);

		// R2D:
		// initialize necessary outputs, convert from Brep to r2d's internal representation...
		r2d_poly polygon;
		r2d_init_poly(&polygon, randverts, nverts); 
		r2d_print(&polygon);
		//... then clip the poly using r2d's function
		r2d_clip(&polygon, planes, nplanes);
		printf("\n");		
		r2d_print(&polygon);
	
		r2d_brep polybrep;
		polybrep.verts = randverts;
		polybrep.nverts = nverts;
		/*
		inplace_clip_poly(&polybrep, planes, nplanes);
		//printf("OUT ");
		info_poly2(&polybrep); */

		r2d_print_brep(&polybrep);
		r2d_clip_brep(&polybrep, planes, nplanes);
		r2d_print_brep(&polybrep);

	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, function: %f\n", elapsed3);
	
}
