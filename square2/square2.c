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
	r2d_real sdists[*nverts], smin, smax, onv;
	r2d_int clipped[MAX_VERTS];
	
	const r2d_real ZERO = 0.0;
	
	
	for(int p = 0; p < nplanes; ++p) {
		onv = *nverts;
		smin = 1.0e30;
		smax = -1.0e30;
		memset(&clipped, 0, sizeof(clipped));

		for(int v = 0; v < onv; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
			if(sdists[v] < smin) { smin = sdists[v]; } 
			if(sdists[v] > smax) { smax = sdists[v]; }
			printf("%f\n", sdists[v]);
			if(sdists[v] < ZERO) { clipped[v] = 1; }
		}

		if(smin >= 0.0) continue;
		if(smax <= 0.0) {
			*nverts = 0;
			return;
		}
		/*for(int v = 0; v < nverts; ++v) {
			if (sdists[v] >= 0) {
				newverts[lastavail] = verts[v];
				lastavail++;
			}
			int vnext = (v+1)%nverts;
			if (sdists[v] * sdists[vnext] < ZERO) {
				r2d_real if1 = sdists[vnext] / (sdists[vnext] - sdists[v]);
				r2d_rvec2 newvert;
				newvert.x = if1*verts[v].x + (1-if1)*verts[vnext].x;
				newvert.y = if1*verts[v].y + (1-if1)*verts[vnext].y;
				newverts[lastavail] = newvert;
				lastavail++;
			}	
		}
		*finalnverts=lastavail;
*/

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

	r2d_plane planes[] = {{{1,0}, -5.0}};
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes

	r2d_rvec2 *randpoly = init_random_poly(&nverts);
	inplace_clip_poly(randpoly, &nverts, planes, nplanes);
	printf("%d\n", nverts);
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
