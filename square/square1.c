#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
//#define REP_TIMES 10000000
//#define REP_TIMES 1
#define REP_TIMES 2

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

// Buffer limit
#define MAX_VERTS 64
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

//takes a poly, clips it, and returns the vertices in the new poly in structure newverts. Incompatible with multiple planes.
void clip_poly(r2d_rvec2* verts, r2d_int nverts, r2d_plane* planes, r2d_int nplanes, r2d_rvec2** newverts, r2d_int* finalnverts) {
	
	//small overallocation for the buffer of the final shape
	//The likely maximum bound is 1.5*nverts; not certain, thus settling for 2*nverts
	*newverts= (r2d_rvec2 *)malloc(MAX_VERTS*sizeof(r2d_rvec2));
	int lastavail = 0;
	r2d_real sdists[nverts];

	const double ZERO = 0.0;
	
	//calculate signed distances
	for(int v = 0; v < nverts; ++v) {
		sdists[v] = planes[0].d + dot(verts[v], planes[0].n);
		//printf("%f\n", sdists[v]);
	}

	//clip the poly and add vertices to newverts buffer
	for(int v = 0; v < nverts; ++v) {
		if (sdists[v] >= 0) {
			(*newverts)[lastavail] = verts[v];
			lastavail++;
		}
		int vnext = (v+1)%nverts;
		if (sdists[v] * sdists[vnext] < ZERO) {
			r2d_real alpha = sdists[vnext] / (sdists[vnext] - sdists[v]);
			r2d_rvec2 newvert;
			newvert.x = alpha*verts[v].x + (1-alpha)*verts[vnext].x;
			newvert.y = alpha*verts[v].y + (1-alpha)*verts[vnext].y;
			(*newverts)[lastavail] = newvert;
			lastavail++;
		}	
	}
	*finalnverts=lastavail;
	//(*newverts) = (r2d_rvec2 *)realloc((*newverts), lastavail*sizeof(r2d_rvec2));
}

//free allocated memory
void free_poly(r2d_rvec2* verts) {
	free(verts);
}

//clips a poly, takes a statically allocated buffer and works in place
void inplace_clip_poly(r2d_rvec2* verts, r2d_int *nverts,
r2d_plane* planes, const r2d_int nplanes) {

	//large static allocation of arrays storing signed distances and kept indices
	r2d_real sdists[MAX_VERTS];
	r2d_int index[MAX_VERTS], newind, oldind, vnext, onv;
	
	const r2d_real ZERO = 0.0;
	
	//loop over each clipping plane
	for(int p = 0; p < nplanes; p++) {
		
		r2d_plane plane=planes[p];
		
		onv = *nverts;
		newind = 0;
		oldind = onv;

		//calculate signed distances
		for(int v = 0; v < onv; ++v) {
			sdists[v] = plane.d + dot(verts[v], plane.n);
		}

		//loop over each vertex
		for(int v = 0; v < onv; ++v) {

			//note the indices of vertices that will be kept
			if (sdists[v] >= 0){
				index[newind++] = v;
			}
			vnext = (v+1)%onv;
			
			//interpolate new vertex if edge crosses
			if (sdists[v] * sdists[vnext] < ZERO) {
				r2d_real alpha = sdists[vnext] / (sdists[vnext] - sdists[v]);
				r2d_rvec2 newvert;
				newvert.x = alpha*verts[v].x + (1-alpha)*verts[vnext].x;
				newvert.y = alpha*verts[v].y + (1-alpha)*verts[vnext].y;
				verts[oldind] = newvert;
				index[newind++] = oldind;
				oldind++;
			}	
		}
		
		//vertices with indices noted in index[] are correct, transfer to buffer
		for(int i = 0; i < newind; ++i) 
			verts[i] = verts[index[i]];
		*nverts = newind;
	}
}

//generate a random poly in a buffer of the right size, for use by clip_poly and r2d_clip
r2d_rvec2* init_random_poly(int *numverts) {

	//generate random number of vertices
	int randint = (rand()%8)+3;
	*numverts = randint;
	
	//allocate the random poly
	r2d_rvec2* verts = (r2d_rvec2 *)malloc(MAX_VERTS * sizeof(r2d_rvec2));

	//randomly generate coordinates and store each set of coordinates as a vertex in the random poly
	for(int i = 0; i < randint; ++i){
		r2d_real randx = RAND_COORD_LIMIT*(double)rand()/(double)RAND_MAX;
		r2d_real randy = RAND_COORD_LIMIT*(double)rand()/(double)RAND_MAX;
		r2d_rvec2 randvert;
		randvert.x = randx;
		randvert.y = randy;
		verts[i] = randvert;	
	}
	return verts;
}

//diagnostic information, printing coordinate information of the random polygon 
void info_poly(r2d_rvec2 *poly, r2d_int nverts) { 
	printf("\n______________________________________________\n");
	
	//number of vertices	
	printf("%d vertices in poly: \n", nverts);
	
	//all vertex coordinates
	printf("Vertex coordinates:\n");
	for (int i = 0; i < nverts; ++i) {
		r2d_rvec2 printvert = poly[i];
		printf("Vertex %d - ", i);
		printf("X: %f ", printvert.x);
		printf("Y: %f\n", printvert.y);
	}
	printf("______________________________________________\n");
}

int main(){

	//number of vertices
	r2d_int nverts;
	srand(time(0));

	// the intersection planes
	r2d_plane planes[] = {{{1,0}, -5.0}};
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes
	// simple timer loop	
	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) 
		{

		// generate a random poly, used by all functions
		r2d_rvec2 *randpoly = init_random_poly(&nverts);
		//info_poly(randpoly, nverts); // information about original random poly

		// OUR RETURNER:
		//initialize necessary outputs, then clip the poly using our prototype function
		r2d_rvec2* finalverts;
		r2d_int finalnverts;
		clip_poly(randpoly, nverts, planes, nplanes, &finalverts, &finalnverts);
		//info_poly(finalverts, finalnverts); // information about returned poly
		//free_poly(finalverts);

		// R2D:
		// initialize necessary outputs, convert from Brep to r2d's internal representation...
		r2d_poly polygon;
		r2d_init_poly(&polygon, randpoly, nverts); 
		
		//... then clip the poly using r2d's function
		r2d_clip(&polygon, planes, nplanes);
		//r2d_print(&polygon);
	
		// OUR INPLACE FUNCTION:
		//clip the poly in place, using our protoype
		info_poly(randpoly, nverts); // information about inplace poly
		inplace_clip_poly(randpoly, &nverts, planes, nplanes);
		info_poly(randpoly, nverts); // information about inplace poly
		
		//free the output memory; we do not need it for now.
		free_poly(randpoly);
	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, function: %f\n", elapsed3);
	
}
