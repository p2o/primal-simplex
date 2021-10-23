Foi criado o algoritmo Primal Simplex como uma classe no arquivo "app/primal_simplex.py"


Para importar a classe faça:<br>
from app.primal_simplex import PrimalSimplex<br>


Para criar uma instância faça:<br>
pl = PrimalSimplex(A,b,c)<br>

onde temos:<br>
A:= é a matriz dos coeficientes das restrições;<br>
b:= é o vetor dos recursos (ou termos independentes);<br>
c:= é o vetor dos custos;<br>


Os métodos da classe são:<br>
solveFi():= é o método que irá encontrar uma base factível para para solução do problema de programação linear;<br>
solve():= é o método que soluciona o problema de programação linear dado;<br>


Os atributos da classe são:<br>
A:= é a matriz dos coeficientes das restrições;<br>
b:= é o vetor dos recursos (ou termos independentes);<br>
c:= é o vetor dos custos;<br>
base:= é o vetor de índices das colunas da matriz A que serão usadas como base;<br>
nbase:= é o vetor de índices das colunas da matriz A que serão usadas como não base;<br>
x:= é o vetor com a solução ótima;<br>
fx:= é o valor ótimo da solução encontrada;<br>

