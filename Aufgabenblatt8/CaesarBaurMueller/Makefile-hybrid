# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-par-hybrid
OBJS = partdiff-par-hybrid.o askparams.o

# Targets ...
all: partdiff-par-hybrid

partdiff-par-hybrid: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-par-hybrid.o: partdiff-par-hybrid.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
