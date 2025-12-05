CC = gcc
CFLAGS = -O3 -march=native -mavx -fopt-info-vec -Wall -DLIKWID_PERFMON
LFLAGS = -lm -llikwid

PROG = cgSolver
MODULES = utils sislin $(PROG)
OBJS = $(addsuffix .o, $(MODULES))

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

cgSolver.o: cgSolver.c sislin.h utils.h
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	@rm -f core *~ *.bak *.o

purge: clean
	@rm -f $(PROG)