Fase 0 - Fundação Conceitual

1. Matriz Simétrica: 
 - Uma matriz A é simétrica se for igual à sua transposta. Ná prática, uma deve ser o espelho da outra em relação à diagonal principal.
 A = |1 3 5| é simétrica mas B = |1 3 5| não é pois B[1][2] != B[2][1].
     |3 9 8|                     |4 9 8|
     |5 8 4|                     |5 8 4|

2. Matriz Positiva Definida
 - Geometricamente, para qualquer vetor de entrada x (exceto o nulo), Ax nunca gira mais que 90 graus.
 - Matematicamente, a operação xT * Ax > 0 deve ser verdadeira para ser Positiva Definida.

3. Matriz Estritamente Diagonalmente Dominante:
 - Em uma linha da matriz A, o módulo do termo da diagonal dominante deve ser maior que a soma dos módulos dos outros elementos da linha.
 - Aplicando para todas o conceito anterior para todas as linhas, se tudo for satisfeito, A é estritamente diagonalmente dominante.

4. Método do Gradiente Conjugado:
 - Objetivo de minimizar uma função, pois encontrar a solução de Ax = b é o mesmo que encontrar o ponto x que satisfaz ∇f(x) = Ax + b = 0.
 - É como o método do gradiente normal, onde você encontra a direção de máxima inclinação negativa e dá um passo nessa direção, repetindo o processo até encontrar um mínimo, porém a partir do segundo cálculo do gradiente, utiliza-se a direção de busca anterior. Por isso o nome conjugado.

