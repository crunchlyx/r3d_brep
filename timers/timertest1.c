
/* 
A simple timing program, giving the time elapsed in ms to run a program a given amount of times defined by the variable REP_TIMES.
*/

#include <time.h>
#include <stdio.h>


//The amount of times the program will be performed. 
#define REP_TIMES 1000000

// A simple program to test the timer. Returns the sum of all elements of a small array.

int addElements(int arry[]){
	int sum;
	int size = sizeof(arry) / sizeof(arry)[0];
	int j = 0;
	for (j; j < size; j++)
		sum+=arry[j];
	return sum; 
}	

//int A[] = {1,2,3,4,5};




//The actual timer, elapsed time given in ms.
int main(){
    clock_t start = clock();
	for (int i = 0; i<REP_TIMES; i++)
		{	
		//code goes here
		}
    clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed);
}
