#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 10000000
//#define REP_TIMES 2
//#define REP_TIMES 1

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

void init_random_poly(r2d_rvec2* verts, r2d_int* nverts) {

	// generate random number of vertices
	//int randint = (rand()%8)+3;

	int randint = 15;
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

void init_nice_square(r2d_rvec2* verts, r2d_int* nverts) {
	*nverts = 4;
	verts[0].x = 1.0;
	verts[0].y = 1.0;
	verts[1].x = 1.0;
	verts[1].y = 9.0;
	verts[2].x = 9.0;
	verts[2].y = 9.0;
	verts[3].x = 9.0;
	verts[3].y = 1.0;
}

int main(){

	r2d_int nverts;
	r2d_int npoly, npoly2;
	srand(time(0));
	r2d_rvec2 randverts[R2D_MAX_VERTS];
	r2d_real momentspos[1];
	r2d_real momentspos2[1];
	r2d_real momentsneg[1];
	r2d_real momentspneg2[1];
	
	// the intersection planes
	r2d_plane planes[] = {{{1,0}, -2.0}, {{1,0}, -5.0},{{1,0}, -8.0}};
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes

	r2d_plane plane = {{1,0}, -5.0};
	// simple timer loop	
	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) {
		
		init_random_poly(randverts, &nverts);

		r2d_brep outpolypos[1];
		r2d_brep outpolyneg[1];
		r2d_brep inpolys[1];
		for (int k = 0; k < 1; ++k){
			for (int j =0; j<nverts; ++j) {
				inpolys[k].verts[j]=randverts[j];
			}
			inpolys[k].nverts = nverts;
		}
		npoly = 1;


		r2d_split_brep(&inpolys, npoly, plane, &outpolypos, &outpolyneg);
		r2d_reduce_brep(&outpolypos, &momentspos2, 1);
		r2d_reduce_brep(&outpolyneg, &momentsneg2, 1);
		/*r2d_brep polybrep;
		for (int j =0; j<nverts; ++j) {
			polybrep.verts[j]=randverts[j];
		}
		polybrep.nverts = nverts;
		r2d_brep polybrep2;
		for (int j =0; j<nverts; ++j) {
			polybrep2.verts[j]=randverts[j];
		}
		polybrep2.nverts = nverts;*/

		// R2D:
		r2d_poly polygon[1];
		r2d_poly r2d_outpos[1];
		r2d_poly r2d_outneg[1];
		r2d_init_poly(&polygon, randverts, nverts); 
		npoly2 = 1;
		
		r2d_split(&polygon, npoly2, plane, &r2d_outpos, &r2d_outneg);
		r2d_reduce(&r2d_outpos, &momentspos, 1);
		r2d_reduce(&r2d_outneg, &momentsneg, 1);

		for(int z = 0; z < 1; ++z){
			if(momentspos[z] != momentspos2[z]){
				printf("component = %d\n", z);
				printf("Moment1 = %f\n", momentspos[z]);
				printf("Moment2 = %f\n", momentspos2[z]);
			}
		}
		for(int z = 0; z < 1; ++z){
			if(momentsneg[z] != momentsneg2[z]){
				printf("component = %d\n", z);
				printf("Moment1 = %f\n", momentsneg[z]);
				printf("Moment2 = %f\n", momentsneg2[z]);
			}
		}

		/*r2d_clip(&polygon, planes, nplanes);
		//r2d_print(&polygon);

		// Our work buffer
		r2d_clip_brep(&polybrep2, planes, nplanes);
		//r2d_print_brep(&polybrep2);*/
	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, function: %f\n", elapsed3);
	
}
