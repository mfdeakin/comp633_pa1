
CC=gcc
CFLAGS = -std=gnu99
LIBS = -lm

debug: CFLAGS += -DMTXDEBUG -g -Wall -pg
debug: LDFLAGS += -pg
debug: pa1

release: CFLAGS += -O3 -floop-interchange -floop-strip-mine -floop-block
release: pa1

pa1: pa1.o matrix.o
	$(CC) $(LDFLAGS) $(LIBS) pa1.o matrix.o -o pa1

clean:
	rm -f *.o pa1
