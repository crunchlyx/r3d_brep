
#include "r3d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define dot(va, vb) (va.x * vb.x + va.y * vb.y + va.z * vb.z)
#define R3D_MAX_DEGREE 4
#define REP_TIMES 100000000

#define R3D_MAX_VERTS_PER_FACE 12
#define MOMENT_TOLERANCE 1.0e-15
#define MOMENT 0
#define MOMENT_ARRAY_SIZE 1

//not proud of this, but we may have to do it
#define R3D_MAX_FACES 16

// mostly for convenience, really
typedef struct {
	r3d_int faceind;
	r3d_int vertind;
} edgeind; 

void printer(r3d_brep* poly) {
	r3d_int* nvertsperface = poly->numvertsperface;
	r3d_int** faceinds = poly->faceinds;
	r3d_rvec3* verts = poly->vertices;
	r3d_int nverts = poly->numvertices;
	r3d_int nfaces = poly->numfaces;

	printf("\n______________________________________________\n");
	printf("%d faces in poly, %d vertices in poly\n", nfaces, nverts);
	for(int f = 0; f < nfaces; ++f) {
		printf("Face %d: ", f);
		for(int i = 0; i < nvertsperface[i]; ++i) {
			printf("%d ", faceinds[f][i]);
		}
		printf("\n");
	}
	printf("\n");
	for(int v = 0; v < nverts; ++v) {
		printf("Vertex %d: %f, %f, %f\n", v, verts[v].x, verts[v].y, verts[v].z);
	}

	printf("______________________________________________\n");
}

void r3d_clip_new(r3d_brep* poly, r3d_brep* newpoly, r3d_plane* planes, r3d_int nplanes) {
	
	r3d_int** const original_faceinds = poly->faceinds;
	const double ZERO = 0.0;
	r3d_real sdists[R3D_MAX_VERTS];

	r3d_int newfaceinds[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE];	
	r3d_int faceinds[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE];	
	r3d_rvec3 newverts[R3D_MAX_VERTS];

	r3d_int nverts = poly->numvertices;
	r3d_rvec3* verts = poly->vertices;
	r3d_int nfaces = poly->numfaces;
	r3d_int* nvertsperface = poly->numvertsperface;

	for(int f = 0; f < nfaces; ++f) {
		for(int i = 0; i < nvertsperface[f]; ++i) {
			faceinds[f][i] = original_faceinds[f][i];
		}
	}

	for(int p = 0; p < nplanes; ++p) {

		r3d_int newnverts = 0;
		r3d_int newnfaces = 0;

		r3d_int newnvertsperface[R3D_MAX_VERTS]={0};

		r3d_int oldv2newv[R3D_MAX_VERTS];
		r3d_int olde2newv[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE]={0};
		r3d_int vnext, nextavail, indcur, indnext, firstnewvertsind, facind, verind, verindnext, crossing_vertex;

		// find the surviving vertices
		for(int v = 0; v < nverts; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
			if (sdists[v] >= ZERO) {
				newverts[newnverts] = verts[v];
				oldv2newv[v] = newnverts;
				newnverts++; 
			}
		}
		
		// prepare for adding crossing vertices
		firstnewvertsind = newnverts;

		// if everything lives on negative side of the clip plane, then bail
		if (newnverts == 0) {
			return;
		}

		// ALLOCATE THE INTEGER MATRIX FOR FINDING EDGE INVERSES
		edgeind invedges[nverts][nverts]; 
		edgeind newv2olde[R3D_MAX_VERTS]; // ???
		edgeind* pe;

		// actually compute the inverse edge matrix
		// the two array indices are the old vertex id's, not the f,i pair
		for(int f = 0; f < nfaces; ++f) {
			for(int v = 0; v < nvertsperface[f]; ++v) {
				indcur = faceinds[f][v];
				indnext = faceinds[f][v == nvertsperface[f]-1 ? 0 : v+1];

				// reverse the indices so this gives the inverse edge, in other
				// words, this f,i edgeinfo is assigned to the opposite edge
				pe = &invedges[indnext][indcur];
				pe->faceind = f;
				pe->vertind = v;
			}
		}
		
		// find the edges that cross and create a new vertex
		for(int f = 0; f < nfaces; ++f) {
			nextavail = 0;
			for(int v = 0; v < nvertsperface[f]; ++v) {
				// check for edge crossings and insert new vertices
				indcur = faceinds[f][v];
				vnext = v==nvertsperface[f]-1 ? 0 : v+1;
				indnext = faceinds[f][vnext];

				//  this edge crosses, and we have not see it before
				if (!(olde2newv[f][v]) && (sdists[indcur] * sdists[indnext]) < ZERO){

						r3d_real alpha = sdists[indnext] / (sdists[indnext] - sdists[indcur]);
						newverts[newnverts].x = alpha*verts[indcur].x + (1-alpha)*verts[indnext].x;
						newverts[newnverts].y = alpha*verts[indcur].y + (1-alpha)*verts[indnext].y;
						newverts[newnverts].z = alpha*verts[indcur].z + (1-alpha)*verts[indnext].z;
					
						// construct edge to new vertex mapping
						olde2newv[f][v]= newnverts;


						r3d_int _f = invedges[indcur][indnext].faceind;
						r3d_int _v = invedges[indcur][indnext].vertind;

						// map the newly created vertex on both this edge and the inverse edge
						// we always use the descending edge
						olde2newv[f][v] = olde2newv[_f][_v] = newnverts;
						if( sdists[indcur]>ZERO){
							// the original edge is descending
							newv2olde[newnverts].faceind = f;
							newv2olde[newnverts].vertind = v;
						} else {
							// the original edge is ascending, 
							newv2olde[newnverts].faceind = _f;
							newv2olde[newnverts].vertind = _v;
						}
					
						// bump the new vertex counter
						newnverts++;				
				}	

				// do the clip 		
				if (sdists[indcur] >= ZERO){
					newfaceinds[newnfaces][nextavail++] = oldv2newv[indcur];
					newnvertsperface[newnfaces]++;
				}	//push oldv2newv[indcur] onto newfaceinds[newnfaces], bump newnvertsperface

				if (olde2newv[f][v] > 0) {
					newfaceinds[newnfaces][nextavail++] = olde2newv[f][v]; 
					newnvertsperface[newnfaces]++;
				}
				// if olde2newv[f][v]>0, that means an edge crossed; push olde2newv[f][v] to newfaceinds[newnfaces]
			}

			if(newnvertsperface[newnfaces] > 0){
				newnfaces++;
			}
		}

		r3d_int marked[R3D_MAX_VERTS] = {0};

		// stitch up the crossing vertices into new faces
		// loop over all new vertices to make sure none are orphaned
		for(int v = firstnewvertsind; v<newnverts; ++v){

			// this begins a new face
			nextavail = 0;
			if (marked[v]) continue;

			// this edge is always descending and crossing by construction
			facind = newv2olde[v].faceind;
			verind = newv2olde[v].vertind;	 
			crossing_vertex = olde2newv[facind][verind];

			do  {
				indcur = faceinds[facind][verind];
				indnext = faceinds[facind][verind == nvertsperface[facind]-1 ? 0 : verind+1];
				// if edge crosses then push
				if (crossing_vertex) {

					// push the new vertex onto the new face
					newfaceinds[newnfaces][nextavail++] = crossing_vertex;
					marked[crossing_vertex] = 1;

					// advance the index in the old face (used to calculate inverse edge)
					verindnext = verind == nvertsperface[facind]-1 ? 0 : verind+1; 	

					// update to the inverse edge
					facind = invedges[indcur][indnext].faceind;
					verind = invedges[indcur][indnext].vertind;

				} // the edge crossed

				// advance the edge (potentially in the new face)
				verind = verind == nvertsperface[facind]-1 ? 0 : verind+1; 			

				// cache the edge to crossing index if it exists
				crossing_vertex = olde2newv[facind][verind];
			
			} while (crossing_vertex!=v);
			newnvertsperface[newnfaces++] = nextavail; 
		}

		// Lets get ready to do it again with another plane; update all old parameters 
		// (may be able to save time on last plane by not doing this step.)
		nverts = newnverts;
		nfaces = newnfaces;
		if(p == nplanes-1){
			break;
		}	
		for(int v = 0; v < newnverts; ++v) {
			verts[v] = newverts[v];
		}	
		for(int f = 0; f < newnfaces; ++f) {
			nvertsperface[f] = newnvertsperface[f];
			for(int i = 0; i < newnvertsperface[f]; ++i) {
				faceinds[f][i] = newfaceinds[f][i];
			}
		}
	}

	// FINAL STEP: We have gone through all planes. Time to pack it up
	
	// malloc a properly sized poly
	r3d_rvec3* finalverts = (r3d_rvec3*) malloc(nverts * sizeof(r3d_rvec3));
	r3d_int* finalnvertsperface = (r3d_int*) malloc(nfaces * sizeof(r3d_int));
	r3d_int** finalfaceinds = (r3d_int**) malloc(nfaces * sizeof(r3d_int*));
	for(int i = 0; i < nfaces; ++i) {
		finalfaceinds[i] = (r3d_int*) malloc(sizeof(r3d_int) * nvertsperface[i]);
	}

	// copy everything
	for(int v = 0; v < nverts; ++v) {
		finalverts[v] = newverts[v];
	}	
	for(int f = 0; f < nfaces; ++f) {
		finalnvertsperface[f] = nvertsperface[f];
		for(int i = 0; i < nvertsperface[f]; ++i) {
			finalfaceinds[f][i] = faceinds[f][i];
		}
	}

	// connect pointers; lets go home
	newpoly->numvertsperface = finalnvertsperface;
	newpoly->vertices = finalverts;
	newpoly->faceinds = finalfaceinds;
	newpoly->numvertices = nverts;
	newpoly->numfaces = nfaces;
}

int main() {
	srand(time(0));

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
	for(int x = 0; x < REP_TIMES; ++x) {
		
		r3d_plane planes[] = {{{0,0,1}, -0.5}, {{0,0,1}, -0.75}};
		r3d_int nplanes = sizeof(planes) / sizeof(planes[0]);

		r3d_poly cube, cubebrep;
		r3d_init_poly(&cube, verts, nverts, faceinds, nvertsperface, nfaces);
		r3d_clip(&cube, planes, nplanes); 

		r3d_brep poly, newpoly;
		poly.numvertices = nverts;
		poly.vertices = verts;
		poly.numvertsperface = nvertsperface;
		poly.numfaces = nfaces;
		poly.faceinds = faceinds;
		r3d_clip_new(&poly, &newpoly, planes, nplanes);

		free(newpoly.vertices);
		for (int f = 0; f < newpoly.numfaces; ++f) {
			free(newpoly.faceinds[f]);
		}
		free(newpoly.faceinds);
		free(newpoly.numvertsperface);
	}
}



