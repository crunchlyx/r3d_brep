
output: square3.o
	gcc -o square3 square3.o -L../../. -lr3d -lm 

square3.o: square3.c
	gcc -c -g  -std=c99 -I ../../. square3.c

clean:
	rm -r *.o core* square3
