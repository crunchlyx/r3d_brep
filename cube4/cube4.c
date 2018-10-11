
#include "r3d.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>

#define dot(va, vb) (va.x * vb.x + va.y * vb.y + va.z * vb.z)
#define R3D_MAX_DEGREE 4
#define REP_TIMES 1

#define R3D_MAX_VERTS_PER_FACE 12
#define MOMENT_TOLERANCE 1.0e-15
#define MOMENT 0
#define MOMENT_ARRAY_SIZE 1

//only for debug, when real version comes around use R3D_MAX_VERTS
#define R3D_MINI_MAX 4

//not proud of this, but we may have to do it
#define R3D_MAX_FACES 16

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

		// if we dont care about preset values, we may toss these out of the loop
		r3d_int newv2newv[R3D_MAX_VERTS];
		r3d_int olde2newv[nverts][nverts];
		memset(olde2newv, 0, sizeof(r3d_int) * nverts * nverts);
		r3d_int oldv2newv[R3D_MAX_VERTS];

		r3d_int nextavail, vcur, vnext, firstnewvertsind, crossing_vertex, dscv, ascv;

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
		
		// find the edges that cross and create a new vertex
		for(int f = 0; f < nfaces; ++f) {
			nextavail = 0;
			
			// a crossing vertex fundamentally can never have index 0, so 0 is fine	
			ascv = 0;
			dscv = 0;
			
			for(int i = 0; i < nvertsperface[f]; ++i) {
				// current vertex
				vcur = faceinds[f][i];

				// if current vertex is positive, push
				if (sdists[vcur] >= ZERO){
					newfaceinds[newnfaces][nextavail++] = oldv2newv[vcur];
					newnvertsperface[newnfaces]++;
				}

				// next vertex
				vnext = faceinds[f][i == nvertsperface[f]-1 ? 0 : i+1];

				bool asc = (sdists[vcur] < ZERO  && sdists[vnext] >= ZERO);
				bool dsc = (sdists[vcur] >= ZERO && sdists[vnext] < ZERO);

				//  this edge crosses and we have not seen it before
				if (asc || dsc) {

					r3d_int crossing_vertex = olde2newv[vnext][vcur];
					
					// if the inverse has not computed the crossing vertex
					if (!crossing_vertex) {

						// compute and create new vertex
						r3d_real alpha = sdists[vnext] / (sdists[vnext] - sdists[vcur]);
						newverts[newnverts].x = alpha*verts[vcur].x + (1-alpha)*verts[vnext].x;
						newverts[newnverts].y = alpha*verts[vcur].y + (1-alpha)*verts[vnext].y;
						newverts[newnverts].z = alpha*verts[vcur].z + (1-alpha)*verts[vnext].z;
						crossing_vertex = newnverts;
						olde2newv[vcur][vnext] = newnverts;
						newnverts++;	
					}

					// push crossing vertex
					newfaceinds[newnfaces][nextavail++] = crossing_vertex;
					newnvertsperface[newnfaces]++;

					// this edge ascends, store crossing vertex for later
					if (asc) {
						ascv = crossing_vertex;
					}
					
					// this edge descends
					else {
						// there are no dangling ascenders, stash vertex
						if (ascv == 0) {
							dscv = crossing_vertex;
						}
						// there is a dangling ascender, assign crossing vertex and clear slate
						else {
							newv2newv[ascv] = crossing_vertex; 
							// clear slate for next ascender
							ascv = 0;
						}
					}
				}	
			}
			
			// there remains one last dangling ascender, assign stashed value
			if (ascv != 0) {
				newv2newv[ascv] = dscv;
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

			vcur = newv2newv[v];

			/*do  {
				// push the new vertex onto the new face
				newfaceinds[newnfaces][nextavail++] = vcur;
				marked[crossing_vertex] = 1;
				vcur = newv2newv[vcur];
			
			} while (vcur!=v);*/
			 while (vcur!=v) {
				// push the new vertex onto the new face
				newfaceinds[newnfaces][nextavail++] = vcur;
				marked[vcur] = 1;
				vcur = newv2newv[vcur];
			}
			newfaceinds[newnfaces][nextavail++] = vcur;
			marked[vcur] = 1;

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

	r3d_int nverts = 4;
	r3d_int nfaces = 4;
	r3d_int nvertsperface[4] = {3,3,3,3};
	r3d_int f0[3] = {0,2,1};
	r3d_int f1[3] = {1,3,0};
	r3d_int f2[4] = {1,2,3};
	r3d_int f3[4] = {3,2,0};
	r3d_int* faceinds[4] = {&f0, &f1, &f2, &f3};
	r3d_rvec3 verts[R3D_MAX_VERTS] = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}; 
	for(int x = 0; x < REP_TIMES; ++x) {
		
		r3d_plane planes[] = {{{0,0,1}, -0.5}};
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

		printer(&newpoly);
		/*free(newpoly.vertices);
		for (int f = 0; f < newpoly.numfaces; ++f) {
			free(newpoly.faceinds[f]);
		}
		free(newpoly.faceinds);
		free(newpoly.numvertsperface);*/
	}
}

