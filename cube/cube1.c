
#include <stdio.h>

// A simple printer that prints all BREP information provided to it. 
void printer(int faces, int* verts, int** faceindices){
	printf("Number of faces: %d\n", faces);
	printf("# of vertices by face: \n");	
	for (int k = 0; k<faces; ++k)
		printf("%d ", *(verts+k));
	printf("\nVertex indices: \n");	
	for (int i = 0; i<faces; ++i){
		for (int j = 0; j<*(verts+i); ++j){
			printf("%d ", *(*(faceindices+i)+j));
		}
		printf("\n");
	}
}

//This is for a cube!
int main(){

	int nfaces = 6; // number of faces in poly. Hardcoded for now.
	int nverts[] = {4,4,4,4,4,4}; // number of vertices per face. Hardcoded for now.
	
	// the imperial face army. holds the brep of each face, hardcoded for now.
	int f0[] = {2,6,5,1};
	int f1[] = {1,5,4,0};
	int f2[] = {4,7,3,0};
	int f3[] = {7,6,2,3};
	int f4[] = {6,7,4,5};
	int f5[] = {3,2,1,0};
	
	int* faceinds[] = {f0,f1,f2,f3,f4,f5}; // the list of all faces as defined by the brep, hardcoded for now.

	printer(nfaces, nverts, faceinds);
}
