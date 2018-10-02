#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 10000000
//#define REP_TIMES 1
//#define REP_TIMES 1

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

// Buffer limit
#define MAX_VERTS 64
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

//free allocated memory
/*void free_poly(r2d_brep* poly) {
	free(poly->verts);
}*/

//index version
void inplace_clip_poly_copy(r2d_brep* poly, r2d_plane* planes, r2d_int nplanes) {

	//large static allocation of arrays storing signed distances and kept indices
	r2d_int index[R2D_MAX_VERTS], newind, oldind, vnext, onv;
	const r2d_real ZERO = 0.0;

	r2d_rvec2 *verts = poly->verts;
	onv = poly->nverts;

	r2d_real sdists[onv];
	
	//loop over each clipping plane
	for(int p = 0; p < nplanes; p++) {
		
		r2d_plane plane=planes[p];
		
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
		//stuff below is valid, don't delete please
		r2d_rvec2 dum[R2D_MAX_VERTS];
		//vertices with indices noted in index[] are correct, transfer to buffer
		for(int i = 0; i < newind; ++i) {
			//printf("%d %d\n",i, index[i]);
			dum[i] = verts[index[i]];
			//printf("%f %f\n", verts[i].x,verts[i].y);
		}

		for(int i = 0; i < newind; ++i) {
			//printf("%d %d\n",i, index[i]);
			verts[i] = dum[i];
			//printf("%f %f\n", verts[i].x,verts[i].y);
		}

		poly->nverts = newind;
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

//diagnostic information, printing coordinate information of the random polygon 
void info_poly(r2d_brep* poly) { 
	printf("\n______________________________________________\n");
	
	r2d_rvec2* verts = poly->verts;
	r2d_int nverts = poly->nverts;
	//number of vertices	
	printf("%d vertices in poly: \n", nverts);
	
	//all vertex coordinates
	printf("Vertex coordinates:\n");
	for (int i = 0; i < nverts; ++i) {
		r2d_rvec2 printvert = verts[i];
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
	r2d_rvec2 randverts[R2D_MAX_VERTS];
	// the intersection planes
	r2d_plane planes[] = {{{1,0}, -5.0}};
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes
	// simple timer loop	
	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) 
		{
		
		// generate a random poly, used by all functions
		init_random_poly(randverts, &nverts);

		// R2D:
		// initialize necessary outputs, convert from Brep to r2d's internal representation...
		/*r2d_poly polygon;
		r2d_init_poly(&polygon, randverts, nverts); 
		//r2d_print(&polygon);
		//... then clip the poly using r2d's function
		r2d_clip(&polygon, planes, nplanes);*/
		//r2d_print(&polygon);
	
		// OUR INPLACE FUNCTION:
		//printf("IN ");

		//clip the poly in place, using our protoype
		r2d_brep polybrep;
		for(int x = 0; x < nverts; ++x) {
			polybrep.verts[x] = randverts[x];
		}
		polybrep.nverts = nverts;
		//info_poly(&polybrep);
		inplace_clip_poly_copy(&polybrep, planes, nplanes);
		//printf("OUT ");
		//info_poly(&polybrep); // information about inplace poly
		
		//free the output memory; we do not need it for now.
		//free_poly(&polybrep);

	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, function: %f\n", elapsed3);
	
}
