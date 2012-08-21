CC=g++
CFLAGS=-c -Wall

P1_SRC=hw4.cpp utility.cpp

P1_OBJ=$(P1_SRC:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

all: app 
		
app: $(P1_OBJ) 
	$(CC) -lgmp -lgmpxx $(P1_OBJ) -o $@

clean:
	rm *.o app 2> /dev/null

