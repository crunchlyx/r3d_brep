#include "r2d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//The amount of times the program will be performed for timer. 
//#define REP_TIMES 10000000
//#define REP_TIMES 2
#define REP_TIMES 1

//Largest random coordinate possibly generated
#define RAND_COORD_LIMIT 10

// Buffer limit
#define MAX_VERTS 64
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)

//free allocated memory
/*void free_poly(r2d_brep* poly) {
	free(poly->verts);
}*/

void r2d_split_brep(r2d_brep* poly, r2d_int npolys, r2d_plane plane, r2d_brep* outpoly_pos, r2d_brep* outpoly_neg){

	// variable declarations
	const r2d_real ZERO = 0.0;	
	r2d_int vnext, nverts, outpos_nverts, outneg_nverts;	
	r2d_rvec2 *verts, *outpos_verts, *outneg_verts;

	for(int p = 0; p < npolys; ++p) {

		nverts = poly[p].nverts;
		verts = poly[p].verts; 
		outpos_nverts=0;
		outneg_nverts=0;
		// direct access to the vertex buffer
		verts = poly[p].verts;
		outpos_verts = outpoly_pos[p].verts;
		outneg_verts = outpoly_neg[p].verts;
		nverts = poly[p].nverts;

		r2d_real sdist_cur = plane.d + dot(verts[0], plane.n);

		// loop over each vertex
		for(int v = 0; v < nverts; ++v) {

			// if on correct side, add vertex coordinates to the next avaiable slot
			// in the buffer
			if (sdist_cur >= 0) {
				outpos_verts[outpos_nverts++] = verts[v];
			}
			if (sdist_cur < 0) {
				 outneg_verts[outneg_nverts++] = verts[v];
			}
		
			vnext = (v+1)%nverts;
			r2d_real sdist_next= plane.d + dot(verts[vnext], plane.n);
	
			// interpolate new vertex if edge crosses and add vertex coordinates to
			// next available slot int the buffer
			if (sdist_cur * sdist_next < ZERO) {
				r2d_real alpha = sdist_next / (sdist_next - sdist_cur);
				outneg_verts[outneg_nverts].x = outpos_verts[outpos_nverts].x = alpha*verts[v].x + (1-alpha)*verts[vnext].x;
				outneg_verts[outneg_nverts].y = outpos_verts[outpos_nverts].y = alpha*verts[v].y + (1-alpha)*verts[vnext].y;
				outpos_nverts++;
				outneg_nverts++;
			}	
			sdist_cur=sdist_next;
		}
		// update value of nverts of the polys
		outpoly_pos[p].nverts = outpos_nverts;
		outpoly_neg[p].nverts = outneg_nverts;
	}
}

void r2d_reduce_brep(r2d_brep* poly, r2d_real* moments, r2d_int polyorder) {

	// var declarations
	r2d_int vcur, vnext, m, i, j, corder;
	r2d_real twoa;
	r2d_rvec2 v0, v1; 

	// direct access to vertex buffer
	r2d_rvec2* vertbuffer = poly->verts; 
	r2d_int* nverts = &poly->nverts; 

	// zero the moments
	for(m = 0; m < R2D_NUM_MOMENTS(polyorder); ++m)
		moments[m] = 0.0;

	if(*nverts <= 0) return;

	// Storage for coefficients
	// keep two layers of the triangle of coefficients
	r2d_int prevlayer = 0;
	r2d_int curlayer = 1;
	r2d_real D[polyorder+1][2];
	r2d_real C[polyorder+1][2];

	// iterate over edges and compute a sum over simplices 
	for(vcur = 0; vcur < *nverts; ++vcur) {

		vnext = (vcur+1)%(*nverts);
		v0.x = vertbuffer[vcur].x;
		v0.y = vertbuffer[vcur].y;
		v1.x = vertbuffer[vnext].x;
		v1.y = vertbuffer[vnext].y;
		twoa = (v0.x*v1.y - v0.y*v1.x); 

		// calculate the moments
		// using the fast recursive method of Koehl (2012)
		// essentially building a set of Pascal's triangles, one layer at a time

		// base case
		D[0][prevlayer] = 1.0;
		C[0][prevlayer] = 1.0;
		moments[0] += 0.5*twoa;

		// build up successive polynomial orders
		for(corder = 1, m = 1; corder <= polyorder; ++corder) {
			for(i = corder; i >= 0; --i, ++m) {
				j = corder - i;
				C[i][curlayer] = 0; 
				D[i][curlayer] = 0;  
				if(i > 0) {
					C[i][curlayer] += v1.x*C[i-1][prevlayer];
					D[i][curlayer] += v0.x*D[i-1][prevlayer]; 
				}
				if(j > 0) {
					C[i][curlayer] += v1.y*C[i][prevlayer];
					D[i][curlayer] += v0.y*D[i][prevlayer]; 
				}
				D[i][curlayer] += C[i][curlayer]; 
				moments[m] += twoa*D[i][curlayer];
			}
			curlayer = 1 - curlayer;
			prevlayer = 1 - prevlayer;
		}
	}

	// reuse C to recursively compute the leading multinomial coefficients
	C[0][prevlayer] = 1.0;
	for(corder = 1, m = 1; corder <= polyorder; ++corder) {
		for(i = corder; i >= 0; --i, ++m) {
			j = corder - i;
			C[i][curlayer] = 0.0; 
			if(i > 0) C[i][curlayer] += C[i-1][prevlayer];
			if(j > 0) C[i][curlayer] += C[i][prevlayer];
			moments[m] /= C[i][curlayer]*(corder+1)*(corder+2);
		}
		curlayer = 1 - curlayer;
		prevlayer = 1 - prevlayer;
	}
}

//generate a random poly in a buffer of the right size, for use by clip_poly and r2d_clip
void init_random_poly(r2d_rvec2* verts, r2d_int* nverts) {

	//generate random number of vertices
	int randint = (rand()%8)+3;

	//int randint = 4;
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

	r2d_int nverts;
	r2d_int npoly;
	srand(time(0));
	r2d_rvec2 randverts[R2D_MAX_VERTS];
	
	// the intersection planes
	//r2d_plane planes[] = {{{1,0}, -5.0}, {{0,1}, -5.0}, {{-1,-1}, -15.0}};
	//int nplanes = sizeof(planes) / sizeof(planes[0]); // number of clipping planes

	r2d_plane plane = {{1,0}, -5.0};
	// simple timer loop	
	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) {
		
		init_random_poly(randverts, &nverts);

		r2d_brep outpolypos[2];
		r2d_brep outpolyneg[2];
		r2d_brep inpolys[2];
		for (int k = 0; k < 2; ++k){
			for (int j =0; j<nverts; ++j) {
				inpolys[k].verts[j]=randverts[j];
			}
			inpolys[k].nverts = nverts;
		}
		npoly = 2;


		r2d_split_brep(&inpolys, npoly, plane, &outpolypos, &outpolyneg);
		for (int l = 0; l < npoly; ++l) {
			r2d_print_brep(&outpolypos[l]);
			r2d_print_brep(&outpolyneg[l]);
		}
		/*
		r2d_brep polybrep;
		for (int j =0; j<nverts; ++j) {
			polybrep.verts[j]=randverts[j];
		}
		polybrep.nverts = nverts;

		r2d_brep polybrep2;
		for (int j =0; j<nverts; ++j) {
			polybrep2.verts[j]=randverts[j];
		}
		polybrep2.nverts = nverts;

		// R2D:
		r2d_poly polygon;
		r2d_init_poly(&polygon, randverts, nverts); 
		//r2d_print(&polygon);
		r2d_clip(&polygon, planes, nplanes);
		r2d_print(&polygon);

		// Our work buffer
		r2d_print_brep(&polybrep2);
		r2d_clip_brep(&polybrep2, planes, nplanes);
		r2d_print_brep(&polybrep2);*/
	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, function: %f\n", elapsed3);
	
}
