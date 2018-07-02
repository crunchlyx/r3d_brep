#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 10000000
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

void clip_poly(r2d_rvec2* verts, r2d_int nverts, r2d_plane* planes, r2d_int nplanes) {
//void clip_poly(r2d_rvec2* verts, r2d_int nverts, r2d_plane* planes, r2d_int nplanes, r2d_rvec2* newverts, r2d_int* nout) {
	
	r2d_real sdists[nverts];
	const double ZERO = 0.0;

	//printf("Plane %d:\n", p); 
	for(int v = 0; v < nverts; ++v) {
		sdists[v] = planes[0].d + dot(verts[v], planes[0].n);
		//printf("%f\n", sdists[v]); 
	}
	//printf("\n");

	//The likely bound is 3/2nverts; not certain.
	r2d_rvec2 newverts[2*nverts];

	int lastavail = 0;
	for(int v = 0; v < nverts; ++v) {
		if (sdists[v] > 0) {
			newverts[lastavail++] = verts[v];
		}
		int vnext = (v+1)%nverts;
		if (sdists[v] * sdists[vnext] < ZERO) {
			r2d_real if1 = sdists[vnext] / (sdists[vnext] - sdists[v]);
			r2d_rvec2 newvert;
			newvert.x = if1*verts[v].x + (1-if1)*verts[vnext].x;
			newvert.y = if1*verts[v].y + (1-if1)*verts[vnext].y;
			newverts[lastavail++] = newvert;
		}	
	}
		/*for(int p = 0; p < sizeof(newverts)/sizeof(newverts)[0]; ++p) {
			printf("Vertex %d - ", p); 			
			printf("X: %f ", newverts[p].x);
			printf("Y: %f\n", newverts[p].y);
		}
		printf("\n");*/
	/*r2d_rvec2 *returnable;
	returnable = (r2d_rvec2 *)malloc(sizeof(newverts));
	if(!returnable) {return NULL;}
	
	for(int i = 0; i < sizeof(newverts)/sizeof(newverts)[0]; ++i)
		returnable[i] = newverts[i];
	
	return returnable;*/
}

r2d_rvec2* init_random_poly(int *numverts) {

	srand(time(0));
	int randint = (rand()%8)+3;
	r2d_rvec2 verts[randint];
	for(int i = 0; i < randint; ++i){
		r2d_real randx = rand()%10;
		r2d_real randy = rand()%10;
		r2d_rvec2 randvert;
		randvert.x = randx;
		randvert.y = randy;
		verts[i] = randvert;	
	}
	*numverts = sizeof(verts)/sizeof(verts)[0];
	r2d_rvec2 *returnable;
	returnable = (r2d_rvec2 *)malloc(randint * sizeof(verts)[0]);
	if(!returnable) {return NULL;}
	
	for(int i = 0; i < randint; ++i) {returnable[i] = verts[i];}
	
	return returnable;
}

int main(){
//the vertices in counterclockwise order	
	/*FOR A SQUARE: {0,0},{1,0},{1,1},{0,1}
	FOR A TRIANGLE: {0,0},{1,0},{0.45,1}
	FOR A HEXAGON: {0,0},{1,0},{2,1},{1,2},{0,2},{-1,2}
	FOR A DODECAGON: {0,0},{1,0},{3,1},{4,3},{4,4},{3,6},{1,7},{0,7},{-2,6},{-3,4},{-3,3},{-2,1}*/
	//FOR A NONCONVEX SHAPE: {0,0},{1,0},{0.45,0.5},{1,1},{0,1}
	
	r2d_int nverts;
	r2d_rvec2 *randpoly = init_random_poly(&nverts);
	r2d_poly polygon;
	r2d_plane planes[] = {{{1,0}, -5.0}}; // the clipping planes
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes

	//diagnostic information, printing all needed information of the random polygon 
	printf("\n______________________________________________\n");
	printf("\nNumber of vertices generated in random poly: %d\n", nverts);
	printf("Randomly generated vertices:\n");
	for (int i = 0; i < nverts; ++i) {
		r2d_rvec2 printvert = randpoly[i];
		printf("X: %f ", printvert.x);
		printf("Y: %f\n", printvert.y);
	}
	printf("\n______________________________________________\n\n");
		



	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) 
		{
			clip_poly(randpoly, nverts, planes, nplanes);
		}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms, us: %f\n", elapsed3); 

    clock_t start = clock();
	for (int p = 0; p<REP_TIMES; ++p)
		{
		r2d_init_poly(&polygon, randpoly, nverts); 
		r2d_clip(&polygon, planes, nplanes);
		}
	clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms, r2d: %f\n", elapsed);
	free(randpoly);
}

