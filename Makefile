    CC     = gcc

    CFLAGS = -O0
    LFLAGS = -lm

     VERIF = verificaEP01
      PROG = labZero
      OBJS = utils.o ZeroFuncao.o

# Lista de arquivos para distribuição
DISTFILES = *.c *.h Makefile polinomio*.dat
DISTDIR = ${USER}

.PHONY: clean purge dist all

%.o: %.c %.h DoubleType.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) $(PROG).o
	$(CC) -o $@ $^ $(LFLAGS)

$(VERIF): $(VERIF).c
	$(CC) -Wno-format -o $@ $<

testeFormato: $(PROG) $(VERIF)
	@cat teste.dat | ./$(PROG) | ./$(VERIF)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG) $(VERIF)

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
