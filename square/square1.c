#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 1000000
//#define REP_TIMES 1000000

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

void clip_poly(r2d_rvec2* verts, r2d_int nverts, r2d_plane* planes, r2d_int nplanes, r2d_rvec2** newverts, r2d_int* finalnverts) {
	
	*newverts= (r2d_rvec2 *)malloc(2*nverts*sizeof(r2d_rvec2));

	//The likely bound is 3/2nverts; not certain.
	int lastavail = 0;
	r2d_real sdists[nverts];
	const double ZERO = 0.0;

	for(int v = 0; v < nverts; ++v) {
		sdists[v] = planes[0].d + dot(verts[v], planes[0].n);
		//printf("%f\n", sdists[v]);
	}

	for(int v = 0; v < nverts; ++v) {
		if (sdists[v] >= 0) {
			(*newverts)[lastavail] = verts[v];
			lastavail++;
		}
		int vnext = (v+1)%nverts;
		if (sdists[v] * sdists[vnext] < ZERO) {
			r2d_real if1 = sdists[vnext] / (sdists[vnext] - sdists[v]);
			r2d_rvec2 newvert;
			newvert.x = if1*verts[v].x + (1-if1)*verts[vnext].x;
			newvert.y = if1*verts[v].y + (1-if1)*verts[vnext].y;
			(*newverts)[lastavail] = newvert;
			lastavail++;
		}	
	}
	*finalnverts=lastavail;
	(*newverts) = (r2d_rvec2 *)realloc((*newverts), lastavail*sizeof(r2d_rvec2));
}

//frees the memory
/*void free_poly(r2d_rvec2 *verts, r2d_int nverts) {
	for(int v=0; v<nverts; v++){
		r2d_rvec2 *freeptr = verts+v;
		free(freeptr);
		printf(" pass");
	}
	printf(" pass");
}*/

r2d_rvec2* init_random_poly(int *numverts) {

	int randint = (rand()%8)+3;
	*numverts = randint;
	r2d_rvec2* verts = (r2d_rvec2 *)malloc(randint * sizeof(r2d_rvec2));

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
	printf("%d vertices in poly: \n", nverts);
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
	for (int i = 0; i<REP_TIMES; ++i) {

		// generate the random poly (used by both)
		r2d_rvec2 *randpoly = init_random_poly(&nverts);

		// print the poly
		//info_poly(randpoly, nverts);

		// clip the poly using our prototype function
		r2d_rvec2* finalverts;
		r2d_int finalnverts;
		clip_poly(randpoly, nverts, planes, nplanes, &finalverts, &finalnverts);

		// print the clipped poly
		//info_poly(finalverts, finalnverts);

		//r2d
		r2d_poly polygon;
		r2d_init_poly(&polygon, randpoly, nverts); 
		r2d_clip(&polygon, planes, nplanes);
		
		//free the poly
		//free_poly(finalverts, finalnverts);
	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, us: %f\n", elapsed3); 
}
