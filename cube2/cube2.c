

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

void random_plane_generator(r3d_plane *empty_plane) {
	
	r3d_int direction = (rand()%3);
	if(direction == 0){
		empty_plane->n.x = 1;
		empty_plane->n.y = 0;
		empty_plane->n.z = 0;
	}
	if(direction == 1){
		empty_plane->n.x = 0;
		empty_plane->n.y = 1;
		empty_plane->n.z = 0;
	}
	if(direction == 2){
		empty_plane->n.x = 0;
		empty_plane->n.y = 0;
		empty_plane->n.z = 1;
	}



	empty_plane->d = (double)rand()/(double)RAND_MAX * (-1);

}
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
void r3d_new_brep(r3d_brep* poly) {

	r3d_rvec3* newverts = (r3d_rvec3*) malloc(R3D_MAX_VERTS * sizeof(r3d_rvec3));
	r3d_int* newnvertsperface = (r3d_int*) malloc(R3D_MAX_VERTS * sizeof(r3d_int));
	r3d_int** newfaceinds = (r3d_int**) malloc(R3D_MAX_VERTS * sizeof(r3d_int*));
	for(int i = 0; i < R3D_MAX_VERTS; ++i) {
		newfaceinds[i] = (r3d_int*) malloc(sizeof(r3d_int) * R3D_MAX_VERTS_PER_FACE);
	}
	r3d_int newnverts = 0;
	r3d_int newnfaces = 0;

	poly->numvertsperface = newnvertsperface;
	poly->vertices = newverts;
	poly->faceinds = newfaceinds;
	poly->numvertices = newnverts;
	poly->numfaces = newnfaces;
}

void r3d_clip_brep(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes) {
//void r3d_clip_brep(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes){
	
	const double ZERO = 0.0;
	r3d_real sdists[R3D_MAX_VERTS];

	for(int p = 0; p < nplanes; ++p) {
		r3d_int nverts = poly->numvertices;
		r3d_rvec3* verts = poly->vertices;
		r3d_int nfaces = poly->numfaces;
		r3d_int* nvertsperface = poly->numvertsperface;
		r3d_int** faceinds = poly->faceinds;

		r3d_int newnverts = 0;
		r3d_int newnfaces = 0;
		r3d_rvec3 newverts[R3D_MAX_VERTS];
		r3d_int newnvertsperface[R3D_MAX_VERTS]={0};
		r3d_int newfaceinds[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE];	

		r3d_int oldv2newv[R3D_MAX_VERTS];
		r3d_int olde2newv[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE]={0};
		r3d_int vnext, nextavail, indcur, indnext, firstnewvertsind, facind, verind, verindnext, crossing_vertex;

		for(int v = 0; v < nverts; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
			if (sdists[v] >= ZERO) {
				newverts[newnverts] = verts[v];
				oldv2newv[v] = newnverts;
				newnverts++; 
			}
		}
		
		firstnewvertsind = newnverts;

		if (newnverts == 0) {
			return;
		}
	
		r3d_int degreecounter[R3D_MAX_VERTS]={0}; // we may only need this to nverts
		edgeind edges[R3D_MAX_VERTS][R3D_MAX_DEGREE];
		edgeind newv2olde[R3D_MAX_VERTS]; // ???
		edgeind i;

		// generate n-skeleton structure
		for(int f = 0; f < nfaces; ++f) {
			for(int v = 0; v < nvertsperface[f]; ++v) {

				// this block creates the one skeleton for each vertex
				indcur = faceinds[f][v];
				edges[indcur][degreecounter[indcur]].faceind = f;
				edges[indcur][degreecounter[indcur]].vertind = v;
				degreecounter[indcur]++;
			}
		}

		// tragically the complete 1-skeleton must be known before we can compute an inverse edge
		// and we need to put a newly created vertex on both an edge and its inverse edge. Because 
		// of this we need to duplicate the loop
		// only need to go to the second to last face (hence the -1, nfaces-1) because by the last face
		// we have seen the all opposite edges already

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

						// calculate the opposite edge on the fly
						r3d_int _f, _v;
						for(int e = 0; e < degreecounter[indnext]; ++e) {
							i = edges[indnext][e];
							if(faceinds[i.faceind][i.vertind==nvertsperface[f]-1?0: i.vertind+1] == indcur ) {
								_f = i.faceind;
								_v = i.vertind;
								break;
							}
						}
						// map the newly created vertex on both this edge and the inverse edge
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

		for(int v = firstnewvertsind; v<newnverts; ++v){

			// this begins a new face
			nextavail = 0;
			if (marked[v]) continue;

			// this edge is always descending and crossing by construction
			facind = newv2olde[v].faceind;
			verind = newv2olde[v].vertind;	 
			crossing_vertex = olde2newv[facind][verind];

			do  {
				
				// if edge crosses then push
				if (crossing_vertex) {

					// push the new vertex onto the new face
					newfaceinds[newnfaces][nextavail++] = crossing_vertex;
					marked[crossing_vertex] = 1;

					// advance the index in the old face (used to calculate inverse edge)
					verindnext = verind == nvertsperface[facind]-1 ? 0 : verind+1; 	

					// find the inverse edge
					for(int e = 0; e < degreecounter[faceinds[facind][verindnext]]; ++e) {
						i = edges[faceinds[facind][verindnext]][e];
						if(faceinds[i.faceind][i.vertind==nvertsperface[facind]-1 ? 0 : i.vertind+1] == faceinds[facind][verind]) {
							facind = i.faceind;
							verind = i.vertind;
							break;
						}
					}
				} // the edge crossed

				// advance the edge (potentially in the new face)
				verind = verind == nvertsperface[facind]-1 ? 0 : verind+1; 			

				// cache the edge to crossing index if it exists
				crossing_vertex = olde2newv[facind][verind];
			
			} while (crossing_vertex!=v);
			newnvertsperface[newnfaces++] = nextavail; 
		}
	}
}

void r3d_clip_exp_inverse(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes) {
//void r3d_clip_brep(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes){
	
	const double ZERO = 0.0;
	r3d_real sdists[R3D_MAX_VERTS];

	for(int p = 0; p < nplanes; ++p) {
		r3d_int const nverts = poly->numvertices;
		r3d_rvec3* const verts = poly->vertices;
		r3d_int const nfaces = poly->numfaces;
		r3d_int* const nvertsperface = poly->numvertsperface;
		r3d_int** const faceinds = poly->faceinds;

		r3d_int newnverts = 0;
		r3d_int newnfaces = 0;
		r3d_rvec3 newverts[R3D_MAX_VERTS];
		r3d_int newnvertsperface[R3D_MAX_VERTS]={0};
		r3d_int newfaceinds[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE];	

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
	}
}
r3d_brep* r3d_clip_tcmalloc(r3d_brep* poly, r3d_brep* newpoly, r3d_plane* planes, r3d_int nplanes) {
//void r3d_clip_brep(r3d_brep* poly, r3d_plane* planes, r3d_int nplanes){
	
	const double ZERO = 0.0;
	r3d_real sdists[R3D_MAX_VERTS];

	for(int p = 0; p < nplanes; ++p) {
		r3d_int nverts = poly->numvertices;
		r3d_rvec3* verts = poly->vertices;
		r3d_int nfaces = poly->numfaces;
		r3d_int* nvertsperface = poly->numvertsperface;
		r3d_int** faceinds = poly->faceinds;

		r3d_new_brep(newpoly);
		r3d_int newnverts = newpoly->numvertices;
		r3d_rvec3* newverts = newpoly->vertices;
		r3d_int newnfaces = newpoly->numfaces;
		r3d_int* newnvertsperface = newpoly->numvertsperface;
		r3d_int** newfaceinds = newpoly->faceinds;
		

		r3d_int oldv2newv[R3D_MAX_VERTS];
		r3d_int olde2newv[R3D_MAX_VERTS][R3D_MAX_VERTS_PER_FACE]={0};
		r3d_int vnext, nextavail, indcur, indnext, firstnewvertsind, facind, verind, verindnext, crossing_vertex;

		for(int v = 0; v < nverts; ++v) {
			sdists[v] = planes[p].d + dot(verts[v], planes[p].n);
			if (sdists[v] >= ZERO) {
				newverts[newnverts] = verts[v];
				oldv2newv[v] = newnverts;
				newnverts++; 
			}
		}
		
		firstnewvertsind = newnverts;

		if (newnverts == 0) {
			return newpoly;
		}
	
		r3d_int degreecounter[R3D_MAX_VERTS]={0}; // we may only need this to nverts
		edgeind edges[R3D_MAX_VERTS][R3D_MAX_DEGREE];
		edgeind newv2olde[R3D_MAX_VERTS]; // ???
		edgeind i;
		edgeind* dum;

		// generate n-skeleton structure
		for(int f = 0; f < nfaces; ++f) {
			for(int v = 0; v < nvertsperface[f]; ++v) {

				// this block creates the one skeleton for each vertex
				indcur = faceinds[f][v];
				dum=&(edges[indcur][degreecounter[indcur]]);
				dum->faceind = f;
				dum->vertind = v;
				degreecounter[indcur]++;
			}
		}

		// tragically the complete 1-skeleton must be known before we can compute an inverse edge
		// and we need to put a newly created vertex on both an edge and its inverse edge. Because 
		// of this we need to duplicate the loop
		// only need to go to the second to last face (hence the -1, nfaces-1) because by the last face
		// we have seen the all opposite edges already

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

						// calculate the opposite edge on the fly
						r3d_int _f, _v;
						for(int e = 0; e < degreecounter[indnext]; ++e) {
							i = edges[indnext][e];
							if(faceinds[i.faceind][i.vertind==nvertsperface[f]-1?0: i.vertind+1] == indcur ) {
								_f = i.faceind;
								_v = i.vertind;
								break;
							}
						}
						// map the newly created vertex on both this edge and the inverse edge
						olde2newv[f][v] = olde2newv[_f][_v] = newnverts;
						if( sdists[indcur] > ZERO){
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

		for(int v = firstnewvertsind; v<newnverts; ++v){

			// this begins a new face
			nextavail = 0;
			if (marked[v]) continue;

			// this edge is always descending and crossing by construction
			facind = newv2olde[v].faceind;
			verind = newv2olde[v].vertind;	 
			crossing_vertex = olde2newv[facind][verind];

			do  {
				
				// if edge crosses then push
				if (crossing_vertex) {

					// push the new vertex onto the new face
					newfaceinds[newnfaces][nextavail++] = crossing_vertex;
					marked[crossing_vertex] = 1;

					// advance the index in the old face (used to calculate inverse edge)
					verindnext = verind == nvertsperface[facind]-1 ? 0 : verind+1; 	

					// find the inverse edge
					for(int e = 0; e < degreecounter[faceinds[facind][verindnext]]; ++e) {
						i = edges[faceinds[facind][verindnext]][e];
						if(faceinds[i.faceind][i.vertind==nvertsperface[facind]-1 ? 0 : i.vertind+1] == faceinds[facind][verind]) {
							facind = i.faceind;
							verind = i.vertind;
							break;
						}
					}
				} // the edge crossed

				// advance the edge (potentially in the new face)
				verind = verind == nvertsperface[facind]-1 ? 0 : verind+1; 			

				// cache the edge to crossing index if it exists
				crossing_vertex = olde2newv[facind][verind];
			
			} while (crossing_vertex!=v);
			newnvertsperface[newnfaces++] = nextavail; 
		}

		// DWS maintains this is redundant
		newpoly->numvertices = newnverts;
		newpoly->vertices = newverts;
		newpoly->numfaces = newnfaces;
		newpoly->numvertsperface = newnvertsperface;
		newpoly->faceinds = newfaceinds;	
	}
	return newpoly;
}

int main() {
	srand(time(0));

	r3d_real askdj;
	r3d_int nverts = 9;
	r3d_int nfaces = 9;
	r3d_int nvertsperface[9] = {4,4,4,4,4,3,3,3,3};
	r3d_int f0[20] = {0,1,5,4};
	r3d_int f1[20] = {1,2,6,5};
	r3d_int f2[20] = {2,3,7,6};
	r3d_int f3[20] = {0,4,7,3};
	r3d_int f4[20] = {0,3,2,1};
	r3d_int f5[20] = {4,5,8};
	r3d_int f6[20] = {5,6,8};
	r3d_int f7[20] = {6,7,8};
	r3d_int f8[20] = {4,8,7};
	r3d_int f9[20];
	r3d_int f10[20];
	r3d_int* faceinds[11] = {&f0,&f1,&f2,&f3,&f4,&f5,&f6,&f7,&f8,&f9,&f10};
	r3d_rvec3 verts[R3D_MAX_VERTS] = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
					 				 {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {0.5,0.5,0.4}};


	for(int x = 0; x < REP_TIMES; ++x) {
		
		askdj = (double)rand()/(double)RAND_MAX * (-.5);
		r3d_plane planes[] = {{{0,0,1}, askdj}};
		r3d_int nplanes = sizeof(planes) / sizeof(planes[0]);

		r3d_poly cube, cubebrep;
		r3d_init_poly(&cube, verts, nverts, faceinds, nvertsperface, nfaces);
		r3d_clip(&cube, planes, nplanes); 

		r3d_brep poly, poly2, newpoly; //r3d_brep poly{verts, nverts, faceinds, nvertsperface, nfaces};
		poly.numvertices = nverts;
		poly.vertices = verts;
		poly.numvertsperface = nvertsperface;
		poly.numfaces = nfaces;
		poly.faceinds = faceinds;
		poly2.numvertices = nverts;
		poly2.vertices = verts;
		poly2.numvertsperface = nvertsperface;
		poly2.numfaces = nfaces;
		poly2.faceinds = faceinds;

		r3d_clip_exp_inverse(&poly, planes, nplanes);

		r3d_clip_tcmalloc(&poly2, &newpoly, planes, nplanes);

		free(newpoly.vertices);
		for (int f = 0; f < R3D_MAX_VERTS; ++f) {
			free(newpoly.faceinds[f]);
		}
		free(newpoly.faceinds);
		free(newpoly.numvertsperface);
	}
}

