CC=g++
CFLAGS=-c -Wall

P1_SRC=hw2_part1.cpp zp_util.cpp gf_util.cpp
P2_SRC=hw2_part2.cpp

P1_OBJ=$(P1_SRC:.cpp=.o)
P2_OBJ=$(P2_SRC:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

all: p1 p2
		
p1: $(P1_OBJ) 
	$(CC) -lgmp -lgmpxx $(P1_OBJ) -o $@

p2: $(P2_OBJ)
	$(CC) -lntl $(P2_OBJ) -o $@

clean:
	rm *.out *.o p1 p2

