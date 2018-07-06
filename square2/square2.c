#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 1000000
//#define REP_TIMES 1000000

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

//The max
#define MAX_VERTS 64

#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

void inplace_clip_poly(r2d_rvec2* verts, r2d_int *nverts, r2d_plane* planes, r2d_int nplanes) {

	//The likely bound is 3/2nverts; not certain.
	r2d_real sdists[MAX_VERTS];
	r2d_int index[MAX_VERTS], newind, oldind, vnext, onv;
	
	const r2d_real ZERO = 0.0;
	
	for(int p = 0; p < nplanes; ++p) {
		onv = *nverts;
		newind = 0;
		oldind = onv;

		for(int v = 0; v < onv; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
		}

		for(int v = 0; v < onv; ++v) {

			if (sdists[v] >= 0){
				index[newind++] = v;
			}
			vnext = (v+1)%onv;

			if (sdists[v] * sdists[vnext] < ZERO) {
				r2d_real if1 = sdists[vnext] / (sdists[vnext] - sdists[v]);
				r2d_rvec2 newvert;
				newvert.x = if1*verts[v].x + (1-if1)*verts[vnext].x;
				newvert.y = if1*verts[v].y + (1-if1)*verts[vnext].y;
				verts[oldind] = newvert;
				index[newind++] = oldind;
				oldind++;
			}	
		}

		for(int i = 0; i < newind; ++i) 
			verts[i] = verts[index[i]];
		*nverts = newind;
	}
}



r2d_rvec2* init_random_poly(int *numverts) {

	int randint = (rand()%8)+3;
	*numverts = randint;
	//r2d_rvec2* verts = (r2d_rvec2 *)malloc(randint * sizeof(r2d_rvec2));
	//static allocation
	r2d_rvec2* verts = (r2d_rvec2 *)malloc(MAX_VERTS * sizeof(r2d_rvec2));

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

//provisional runner
int main() {
	r2d_int nverts;
	srand(time(0));

	/*
	FOR A SQUARE: {0,0},{1,0},{1,1},{0,1}
	FOR A TRIANGLE: {0,0},{1,0},{0.45,1}
	FOR A HEXAGON: {0,0},{1,0},{2,1},{1,2},{0,2},{-1,2}
	FOR A DODECAGON: {0,0},{1,0},{3,1},{4,3},{4,4},{3,6},{1,7},{0,7},{-2,6},{-3,4},{-3,3},{-2,1}
	FOR A NONCONVEX SHAPE: {0,0},{1,0},{0.45,0.5},{1,1},{0,1}*/
	r2d_plane planes[] = {{{1,0}, -5.0}, {{0,1}, -5.0}};
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes

	r2d_rvec2 *randpoly = init_random_poly(&nverts);
	info_poly(randpoly, nverts);
	inplace_clip_poly(randpoly, &nverts, planes, nplanes);
	info_poly(randpoly, nverts);

}


/*int main(){

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
		r2d_poly polygon;
		r2d_init_poly(&polygon, randpoly, nverts); 
		r2d_clip(&polygon, planes, nplanes);


		// clip the poly using our prototype function
		inplace_clip_poly(randpoly, nverts, planes, nplanes);

		// print the clipped poly
		//info_poly(finalverts, finalnverts);

		//r2d


		
		//free the poly
		//free_poly(finalverts, finalnverts);
	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, us: %f\n", elapsed3); 
}*/
