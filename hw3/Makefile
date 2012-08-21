CC=g++
CFLAGS=-c -Wall

P1_SRC=main.cpp

P1_OBJ=$(P1_SRC:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

all: p1 
		
p1: $(P1_OBJ) 
	$(CC) -lgmp -lgmpxx $(P1_OBJ) -o $@

clean:
	rm *.out *.o p1 

