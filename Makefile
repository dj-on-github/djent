CC = gcc
CFLAGS = -I/usr/local/include -m64 -g -Wall 
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm

djent: djent.o chisq.o
	$(CC) $(CFLAGS) $(LDFLAGS) chisq.o djent.o -o djent $(LDLIBS)

djent.o: djent.c
	$(CC) -c $(CFLAGS) -o djent.o djent.c

install:
	cp djent /usr/local/bin

clean:
	rm djent.o
	rm djent

