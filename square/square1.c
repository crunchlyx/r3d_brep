#include "r2d.h"
#include <stdio.h>
#include <time.h>

//The amount of times the program will be performed for timer. 
#define REP_TIMES 100000000
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

/*typedef struct r2d_brep
{
	r2d_rvec2 *verts;
	r2d_int nverts;
	
} r2d_brep;*/
//This is for a square!
/*void r2d_init_brep(r2d_poly *poly, r2d_brep brep*) {
	r2d_int nverts = poly->nverts;
	r2d_vertex verts = poly->verts;
	r2d
} */

int main(){

	r2d_rvec2 verts[] = {{0,0},{1,0},{1,1},{0,1}}; //the vertices in counterclockwise order
	r2d_int nverts = sizeof(verts)/sizeof(verts)[0]; //number of verts

	r2d_poly square; // the square
	r2d_plane planes[] = {{{1,0}, -0.5}}; // the clipping planes
	int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes
	r2d_int vertorder[] = {0,1,2,3};
	
	//timer
    /*clock_t start = clock();
	for (int i = 0; i<REP_TIMES; i++)
		{

		r2d_init_poly(&square, verts, nverts); 
		r2d_clip(&square, planes, nplanes);
		}
	clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed);*/

	/*r2d_print(&square);
	printf("\n");
	r2d_clip(&square, planes, nplanes);
	r2d_print(&square);*/

	r2d_real sdists[nverts];
	

	const double ZERO = 0.0;

	for(int p = 0; p < nplanes; ++p) {
		printf("Plane %d:\n", p); 
		for(int v = 0; v < nverts; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
			printf("%f\n", sdists[v]); 
		}
		printf("\n");
	}

	clock_t start2 = clock();
	for (int i = 0; i<REP_TIMES; i++) {
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
			printf("X: %f ", newverts[p].x);
			printf("Y: %f\n", newverts[p].y);
		}*/
	} 
	clock_t stop2 = clock();
	double elapsed2 = (double)(stop2 - start2) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed2); 

    clock_t start = clock();
	for (int i = 0; i<REP_TIMES; i++)
		{
		r2d_init_poly(&square, verts, nverts); 
		r2d_clip(&square, planes, nplanes);
		}
	clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed);
}

