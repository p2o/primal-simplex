from app.primal_simplex import PrimalSimplex

'''
# nao factivel
A = [[1,1,1,0],[2,3,0,-1]]
b = [4,18]
c = [-3,4,0,0]

A = [[-8,3,1,0,0,0],[3,5,0,-1,0,0],[3,-4,0,0,-1,0],[1,0,0,0,0,1]]
b = [24,21,6,3]
c = [-2,-2,0,0,0,0]


# regiao factivel infinita
A = [[-1,1,1,0],[-1/2,1,0,1]]
b = [1,2]
c = [-2,-2,0,0]


# multiplas solucoes otimas
A = [[1,13/2,-1,0,0],[2,1,0,-1,0],[5,4,0,0,1]]
b = [5,4,20]
c = [2,1,0,0,0]


# solucao ótima
A = [[-3,4,1,0,0],[1,-1,0,1,0],[1,1,0,0,1]]
b = [12,4,6]
c = [-1,-3,0,0,0]

A = [[1,1,1,0,0,0,0],[6,2,0,-1,0,0,0],[1,5,0,0,-1,0,0],[1,0,0,0,0,1,0],[0,1,0,0,0,0,1]]
b = [4,8,4,3,3]
c = [2,3,0,0,0,0,0]
'''

# solucao degenerada
A = [[1,1,1,0,0,0],[1,2,0,1,0,0],[2,1,0,0,1,0],[1,-2,0,0,0,1]]
b = [10,15,15,1]
c = [-3,-5,0,0,0,0]



pl = PrimalSimplex(A,b,c)

print("Matriz A:\n",pl.A)
print("Vetor b:\n",pl.b)
print("Vetor c:\n",pl.c)


pl.solveFi()

print("Base:",pl.base)
print("Não base:",pl.nbase)


pl.solve()

print("Vetor x:",pl.x)
print("Valor ótimo:",pl.fx)
