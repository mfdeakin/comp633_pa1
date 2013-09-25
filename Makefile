
CC=gcc
CFLAGS = -std=gnu99
LIBS = -lm -lsqlite3

debug: CFLAGS += -DMTXDEBUG -Wall -g
debug: pa1

release: CFLAGS += -O3 -floop-strip-mine -floop-block -floop-interchange -msse -msse3 -mmmx
release: pa1

pa1: pa1.o matrix.o
	$(CC) $(LDFLAGS) $(LIBS) pa1.o matrix.o -o pa1

clean:
	rm -f *.o pa1
