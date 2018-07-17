
#include "r3d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define dot(va, vb) (va.x * vb.x + va.y * vb.y + va.z * vb.z)

void r3d_clip_brep(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes){
	
	const double ZERO = 0.0;
	r3d_int nverts = poly->numvertices;
	r3d_rvec3* verts = poly->vertices;
	r3d_int nfaces = poly->numfaces;
	r3d_int* nvertsperface = poly->numvertsperface;
	r3d_int** faceinds = poly->faceinds;

	r3d_real sdists[R3D_MAX_VERTS];

	for(int v = 0; v < nverts; ++v) {
		sdists[v] = planes[0].d + dot(verts[v], planes[0].n);
	}
	
	r3d_int vnext;
	for(int f = 0; f < nfaces; ++f) {
		for(int v = 0; v < nvertsperface[f]; ++v) {
			vnext = (v+1)%nvertsperface[f];
			if (sdists[faceinds[f][v]] * sdists[faceinds[f][vnext]] < ZERO){
				printf("Cross @ face %d\n", f);
			}
		}
	}		


}

int main() {

	r3d_int nverts = 8;
	r3d_int nfaces = 6;
	r3d_int nvertsperface[6] = {4,4,4,4,4,4};
	r3d_int f0[4] = {0,1,5,4};
	r3d_int f1[4] = {1,2,6,5};
	r3d_int f2[4] = {2,3,7,6};
	r3d_int f3[4] = {0,4,7,3};
	r3d_int f4[4] = {0,3,2,1};
	r3d_int f5[4] = {4,5,6,7};
	r3d_int* faceinds[6] = {&f0, &f1, &f2, &f3, &f4, &f5};
	r3d_rvec3 verts[R3D_MAX_VERTS] = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
					 				 {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}}; 

	r3d_plane planes[] = {{{1,0,0}, -0.5}};
	r3d_int nplanes = sizeof(planes) / sizeof(planes[0]);

	r3d_poly cube;
	//r3d_init_poly(&cube, verts, nverts, faceinds, nvertsperface, nfaces);
	//r3d_clip(&cube, planes, nplanes); 
	//r3d_print(&cube);

	r3d_brep poly;
	poly.numvertices = nverts;
	poly.vertices = verts;
	poly.faceinds = faceinds;
	poly.numvertsperface = nvertsperface;
	poly.numfaces = nfaces;

	r3d_clip_brep(&poly, planes, nplanes);


}
