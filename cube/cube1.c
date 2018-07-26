
#include "r3d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define dot(va, vb) (va.x * vb.x + va.y * vb.y + va.z * vb.z)
#define R3D_MAX_DEGREE 3
#define REP_TIMES 1

// mostly for convenience, really
typedef struct {
	r3d_int faceind;
	r3d_int vertind;
} edgeind;

void r3d_clip_brep(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes){
	
	const double ZERO = 0.0;
	r3d_real sdists[R3D_MAX_VERTS];

	for(int p = 0; p < nplanes; ++p) {
		r3d_int nverts = poly->numvertices;
		r3d_rvec3* verts = poly->vertices;
		r3d_int nfaces = poly->numfaces;
		r3d_int* nvertsperface = poly->numvertsperface;
		r3d_int** faceinds = poly->faceinds;

		r3d_int newnvertsperface[R3D_MAX_VERTS];
		r3d_rvec3 newverts[R3D_MAX_VERTS];
		r3d_int newfaceinds[R3D_MAX_VERTS][R3D_MAX_VERTS];
		r3d_int newnverts = 0;
		r3d_int newnfaces = 0;

		r3d_int newvertind = 0;
		
		r3d_int mapper[R3D_MAX_VERTS], vertindex[R3D_MAX_VERTS], vnext, nextavail, lastvertind, facindcur, facindnext;

		for(int v = 0; v < nverts; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
			if (sdists[v] >= ZERO) {
				newverts[newnverts++] = verts[v];
				vertindex[newvertind++] = v; 
			}
		}
	
		
		r3d_int degreecounter[R3D_MAX_VERTS]={0}; // we may only need this to nverts
		edgeind edges[R3D_MAX_VERTS][R3D_MAX_DEGREE];
		edgeind inverse_edges[nfaces][R3D_MAX_VERTS]; // find something better than MAX_VERTS

		// generate n-skeleton structure
		for(int f = 0; f < nfaces; ++f) {
			for(int v = 0; v < nvertsperface[f]; ++v) {
				r3d_int indcur = faceinds[f][v];
				edges[indcur][degreecounter[indcur]].faceind = f;
				edges[indcur][degreecounter[indcur]].vertind = v;
				degreecounter[indcur]++;
				printf("%d", faceinds[f][v]);
			}
			printf("\n");
		}

		// generate inverse edge structure
		for(int f = 0; f < nfaces; ++f) { // inverse edge calculation
			for(int v = 0; v < nvertsperface[f]; ++v) {
				vnext = (v+1)%nvertsperface[f];
				r3d_int nextind = faceinds[f][vnext];
				for(int e = 0; e < R3D_MAX_DEGREE; ++e) {
					edgeind nextedge = edges[nextind][e];
					if(faceinds[nextedge.faceind][(((nextedge.vertind)+1)%nvertsperface[f])] == (faceinds[f][v]) ) {
						inverse_edges[f][v].faceind = nextedge.faceind;
						inverse_edges[f][v].vertind = nextedge.vertind;
					}
				}
			}
		}

		// we can probably avoid doing 128x128
		r3d_int crossed[R3D_MAX_VERTS][R3D_MAX_VERTS] = {0}; // for marking edges that we have already crossed
		r3d_int insertverts[R3D_MAX_VERTS][R3D_MAX_VERTS]; // for storing vertices previously calculated
		lastvertind = nverts;
		r3d_int survface = 0;

		for(int f = 0; f < nfaces; ++f) {
		nextavail = 0;
			for(int v = 0; v < nvertsperface[f]; ++v) {
				vnext = (v+1)%nvertsperface[f];
				facindcur = faceinds[f][v];
				facindnext = faceinds[f][vnext];
				
				//if (crossed[f][v] ==  1) continue;
				if (sdists[facindcur] >= ZERO) {
					newfaceinds[survface][nextavail++] = facindcur;
					newnvertsperface[survface]++;
				}

				if (sdists[facindcur] * sdists[facindnext] < ZERO && crossed[f][v] == 0) {

					newnvertsperface[survface]++;
					if ( crossed[inverse_edges[f][v].faceind][inverse_edges[f][v].vertind] == 1 ) {
						
						newfaceinds[survface][nextavail++] = insertverts[f][v];
						crossed[f][v] = 1;
					}
					if ( crossed[inverse_edges[f][v].faceind][inverse_edges[f][v].vertind] == 0 ) {
						// interpolation
						r3d_real alpha = sdists[facindnext] / (sdists[facindnext] - sdists[facindcur]);
						verts[lastvertind].x = alpha*verts[facindcur].x + (1-alpha)*verts[facindnext].x;
						verts[lastvertind].y = alpha*verts[facindcur].y + (1-alpha)*verts[facindnext].y;
						verts[lastvertind].z = alpha*verts[facindcur].z + (1-alpha)*verts[facindnext].z;
						newverts[newnverts++] = verts[lastvertind]; 
						vertindex[newvertind++] = lastvertind;

						newfaceinds[survface][nextavail++] = lastvertind;
						crossed[f][v] = 1;
						insertverts[inverse_edges[f][v].faceind][inverse_edges[f][v].vertind] = lastvertind;
						lastvertind++;
					}
				}
			}
			if (newnvertsperface[survface] != 0) { 
				survface++;
			}
		}

		for (int v; v < newnverts; ++v) {
			mapper[vertindex[v]] = v;
		}

		for (int f = 0; f < survface; ++f) {
			for (int v = 0; v < newnvertsperface[f]; ++v) {
				newfaceinds[f][v] = mapper[newfaceinds[f][v]];
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

	clock_t start3 = clock();
	for (int i = 0; i<REP_TIMES; ++i) {
		r3d_clip_brep(&poly, planes, nplanes);
	}
	clock_t stop3 = clock();
	double elapsed3 = (double)(stop3 - start3) * 1000 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms, function: %f\n", elapsed3);
}

