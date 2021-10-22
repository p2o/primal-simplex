Foi criado o algoritmo Primal Simplex como uma classe no arquivo "app/primal_simplex.py"


Para importar a classe faça:
from app.primal_simplex import PrimalSimplex


Para criar uma instância faça:
pl = PrimalSimplex(A,b,c)

onde temos:
A:= é a matriz dos coeficientes das restrições;
b:= é o vetor dos recursos (ou termos independentes);
c:= é o vetor dos custos;


Os métodos da classe são:
solveFi():= é o método que irá encontrar uma base factível para para solução do problema de programação linear;
solve():= é o método que soluciona o problema de programação linear dado;


Os atributos da classe são:
A:= é a matriz dos coeficientes das restrições;
b:= é o vetor dos recursos (ou termos independentes);
c:= é o vetor dos custos;
base:= é o vetor de índices das colunas da matriz A que serão usadas como base
nbase:= é o vetor de índices das colunas da matriz A que serão usadas como não base
x:= é o vetor com a solução ótima
fx:= é o valor ótimo da solução encontrada

