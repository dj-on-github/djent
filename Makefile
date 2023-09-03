CC = gcc
CFLAGS = -I/usr/local/include -m64 -g -Wall
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm -lgmp -lmpfr

all: djent

mathy_things.o: mathy_things.c mathy_things.h
	$(CC) -c $(CFLAGS) -o mathy_things.o mathy_things.c

filename_parse.o: filename_parse.c filename_parse.h
	$(CC) -c $(CFLAGS) -o filename_parse.o filename_parse.c

markov2p.o: markov2p.c markov2p.h
	$(CC) -c $(CFLAGS) -o markov2p.o markov2p.c

djent.o: djent.c markov2p.h filename_parse.h mathy_things.h
	$(CC) -c $(CFLAGS) -o djent.o djent.c

djent: djent.o markov2p.o filename_parse.o mathy_things.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) mathy_things.o filename_parse.o markov2p.o djent.o -o djent 

install:
	cp djent /usr/local/bin

clean:
	rm -f markov2p.o
	rm -f djent.o
	rm -f djent

