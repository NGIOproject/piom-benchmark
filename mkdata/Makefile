CC=cc
CFLAGS=-g

DEPS=Makefile arralloc.h

OBJ=mkinfile.o arralloc.o

all: mkinfile
.PHONY:all

%.o: %.c $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS)

mkinfile: $(OBJ) 
	$(CC) -o $@ $(OBJ) $(CFLAGS)

clean:
	rm -fr *.o mkinfile
