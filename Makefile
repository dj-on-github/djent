CC = gcc
CFLAGS = -I/usr/local/include -m64 -g -Wall
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm -lgmp -lmpfr

all: djent

markov2p.o: markov2p.c markov2p.h
	$(CC) -c $(CFLAGS) -o markov2p.o markov2p.c

djent.o: djent.c
	$(CC) -c $(CFLAGS) -o djent.o djent.c

djent: djent.o markov2p.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS)  markov2p.o djent.o -o djent 

install:
	cp djent /usr/local/bin

clean:
	rm -f markov2p.o
	rm -f djent.o
	rm -f djent

