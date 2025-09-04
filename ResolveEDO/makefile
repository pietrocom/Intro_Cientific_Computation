# Makefile para o projeto de EDO

# Compilador
CC = gcc

# Flags do compilador: -O2 para otimização, -Wall para mostrar todos os avisos, -std=c99 para compatibilidade
CFLAGS = -O2 -Wall -std=c99 -D_POSIX_C_SOURCE=199309L

# Flags do linker: -lm para a biblioteca matemática (sqrt, exp), -lrt para a biblioteca de tempo real (clock_gettime)
LDFLAGS = -lm -lrt

# Nome do executável final
TARGET = edo_solver

# Lista de arquivos fonte (.c)
SRCS = edo.c gaussSeidel_EqDiff.c utils.c

# Gera a lista de arquivos objeto (.o) a partir da lista de fontes
OBJS = $(SRCS:.c=.o)

# --- Regras ---

# Regra padrão (executada ao rodar "make")
# Declara que o objetivo principal é construir o TARGET
.PHONY: all
all: $(TARGET)

# Regra para linkar os arquivos objeto e criar o executável final
# $@ representa o nome do alvo (edo_solver)
# $^ representa a lista de todas as dependências (os arquivos .o)
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Regra genérica para compilar um arquivo .c em um arquivo .o
# -c diz ao compilador para apenas compilar, sem linkar
# $< representa a primeira dependência (o arquivo .c)
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Regra para limpar os arquivos gerados
# rm -f remove os arquivos sem pedir confirmação
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)