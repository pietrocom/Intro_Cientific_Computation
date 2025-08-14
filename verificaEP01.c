#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MSG_ERRO "\n==> Formato de saída INCORRETO !\n\n"
#define MSG_OK "\n==> Formato de saída CORRETO !\n\n"

int main() {

  int n, k, it;
  double raiz, valorcrit, tempo;
  double t;
  char label[10];

  while ( (k=scanf("%s\n", label)) != EOF) {
    if (k != 1) {
      printf(MSG_ERRO);
      exit(1);
    }
  
    printf("%s\n", label);

    for (int i=0; i < 3; ++i) {
      if ( scanf("%s %lf %lf %d %lf \n", label, &raiz, &valorcrit, &it, &tempo) != 5 ){
	printf(MSG_ERRO);
	exit(1);
      }
      
      if ( strcmp("bissec", label) ){
	printf(MSG_ERRO);
	exit(2);
      }

      printf("%s %+.15e %.15e %3d %.8e \n", label, raiz, valorcrit, it, tempo);
    }
    
    for (int i=0; i < 3; ++i) {
      if ( scanf("%s %lf %lf %d %lf \n", label, &raiz, &valorcrit, &it, &tempo) != 5 ){
	printf(MSG_ERRO);
	exit(3);
      }
      
      if ( strcmp("newton", label) ){
	printf(MSG_ERRO);
	exit(4);
      }

      printf("%s %+.15e %.15e %3d %.8e \n", label, raiz, valorcrit, it, tempo);
    }

    printf(MSG_OK);
  
  }
}
