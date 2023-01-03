CC = mpicc
CFLAGS = -std=c11 -pedantic -Wall -Wextra -Og -ggdb -gdwarf-4 -g
LFLAGS = $(CFLAGS)
TGTS = timempi

all: $(TGTS)

$(TGTS): %: %.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(TGTS)
