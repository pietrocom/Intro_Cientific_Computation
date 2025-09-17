# PROGRAMA
PROG = resolveEDO
OBJS = $(PROG).o edo.o utils.o
VERIF = verificaEP02

# Compilador
CC = gcc

# Opções para compilação com LIKWID. 
# ATENÇÃO: pode ser necessário adicionar os caminhos para os includes e libs do LIKWID,
# por exemplo: -I/usr/local/include -L/usr/local/lib
CFLAGS = -O0 -DLIKWID_PERFMON -I/usr/include
LFLAGS = -lm -llikwid

# Lista de arquivos para distribuição.
DISTFILES = *.c *.h LEIAME Makefile *.dat
DISTDIR = pc24

.PHONY: clean purge dist all

all: $(PROG)

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(VERIF): $(VERIF).c
	$(CC) -Wno-format -o $@ $<

testeFormato: $(PROG) $(VERIF)
	@echo "--- Testando o formato da saida com 'teste.dat' ---"
	@cat teste.dat | ./$(PROG) | ./$(VERIF)

clean:
	@echo "Limpando arquivos temporários ..."
	@rm -f *~ *.bak *.o

purge:  clean
	@echo "Limpando todos os arquivos gerados ..."
	@rm -f $(PROG) $(VERIF)
	@rm -f $(DISTDIR).tgz

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@mkdir -p $(DISTDIR)
	@cp -p $(DISTFILES) $(DISTDIR)/
	@tar -czvf $(DISTDIR).tgz $(DISTDIR)
	@rm -rf $(DISTDIR)