CC = gcc

CFLAGS = -O0 -Wall 
LFLAGS = -lm

PROG = cgSolver

MODULES = utils \
          sislin \
          $(PROG)

OBJS = $(addsuffix .o, $(MODULES))
SRCS = $(addsuffix .c, $(MODULES)) $(addsuffix .h, $(MODULES))

all: $(PROG)

# Regra de linkagem para criar o programa final a partir dos objetos
$(PROG): $(OBJS)
	@echo "==> Linkando o executável: $@"
	$(CC) -o $@ $^ $(LFLAGS)

# Regra para compilar em modo de depuração.
debug: CFLAGS += -g -D__DEBUG__
debug: $(PROG)

cgSolver.o: cgSolver.c sislin.h utils.h
	@echo "==> Compilando: $<"
	$(CC) -c $(CFLAGS) $< -o $@

sislin.o: sislin.c sislin.h utils.h
	@echo "==> Compilando: $<"
	$(CC) -c $(CFLAGS) $< -o $@

utils.o: utils.c utils.h
	@echo "==> Compilando: $<"
	$(CC) -c $(CFLAGS) $< -o $@

DISTFILES = *.c *.h Makefile LEIAME
DISTDIR = pc24

.PHONY: clean purge dist all

clean:
	@echo "Limpando sujeira ....."
	@rm -f core *~ *.bak *.o

purge: clean
	@echo "Fazendo a faxina ....."
	@rm -f $(PROG)

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@mkdir -p $(DISTDIR)
	@cp $(DISTFILES) $(DISTDIR)/
	@tar -czvf $(DISTDIR).tgz $(DISTDIR)
	@rm -rf $(DISTDIR)