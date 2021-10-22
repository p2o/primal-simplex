import numpy as np

def simplexFii(matrizA,bases,naoBase,vetorC,vetorB):
    B = matrizA[:, bases]
    N = matrizA[:, naoBase]
    cb = vetorC[bases]
    cn = vetorC[naoBase]
    b = vetorB

    x = np.linalg.solve(B,b)
    lamb = np.linalg.solve(np.transpose(B),cb)

    hat_c = []
    for i in range(N.shape[1]):
        hat_c.append(cn[i,0] - np.dot(np.transpose(lamb),N[:,i])[0])

    hat_c_min_i = np.argmin(hat_c)
    entra_na_base = naoBase[hat_c_min_i]

    if hat_c[hat_c_min_i] >= 0 and all(x>=0):
        x_r = [0] * matrizA.shape[1]
        for i,j in enumerate(bases):
            x_r[j] = x[i][0]
        return [1,x_r] #solução ótima

    y = np.linalg.solve(B,N[:,hat_c_min_i])

    if max(y) <= 0:
        return [2,x] #sem solução ótima finita

    hat_epi = []
    hat_epi_i = []
    for i in range(len(y)):
        if y[i] > 0:
            hat_epi.append(x[i,0] / y[i])
            hat_epi_i.append(i)

    hat_epi_min_i = hat_epi_i[np.argmin(hat_epi)]
    sai_da_base = bases[hat_epi_min_i]

    bases[hat_epi_min_i] = entra_na_base
    naoBase[hat_c_min_i] = sai_da_base

    return [0,bases,naoBase]

class PrimalSimplex:

    def __init__(self,A,b,c):
        self.A = np.asarray(A)
        self.b = np.asarray(b).reshape(len(b),1)
        self.c = np.asarray(c).reshape(len(c),1)

    def solveFi(self):
        tamanho_base = self.A.shape[0]
        self.base = [-1] * tamanho_base
        self.nbase = []
        matriz_B_to_find = []
        base_i = 0
        id_columns = []
        for i in range(tamanho_base):
            matriz_B_to_find.append([0] * tamanho_base)
            matriz_B_to_find[i][i] = 1
        for col in range(self.A.shape[1]):
            for n_base in range(len(matriz_B_to_find)):
                if all(self.A[:, col] == matriz_B_to_find[n_base]):
                    self.base[base_i] = col
                    id_columns.append(n_base)
                    base_i += 1
                    break
            if base_i >= tamanho_base:
                self.nbase = [i for i in range(self.A.shape[1]) if i not in self.base]
                break

        if base_i < tamanho_base:
            columns_to_add = np.asarray([matriz_B_to_find[i] for i in range(tamanho_base) if i not in id_columns]).T
            A_aux = np.concatenate((self.A,columns_to_add), axis=1)
            c_aux = np.asarray([0] * self.A.shape[1] + [1] * columns_to_add.shape[1]).reshape(A_aux.shape[1],1)
            while base_i < tamanho_base:
                self.base[base_i] = A_aux.shape[1] - tamanho_base + base_i
                base_i += 1
            self.nbase = [i for i in range(A_aux.shape[1]) if i not in self.base]
            x = [0]
            while x[0] == 0:
                x = simplexFii(A_aux,self.base,self.nbase,c_aux,self.b)
                if x[0] == 0:
                    self.base = x[1]
                    self.nbase = x[2]
            for i in range(self.A.shape[1],A_aux.shape[1]):
                if i in self.nbase:
                    self.nbase.remove(i)

        return [self.base, self.nbase]

    def solve(self):
        if not hasattr(self, "base"):
            self.solveFi()
        if self.A.shape[1] < len(self.base) + len(self.nbase):
            print("Problema infactivel")
        else:
            x = [0]
            while x[0] == 0:
                x = simplexFii(self.A,self.base,self.nbase,self.c,self.b)
                if x[0] == 0:
                    self.base = x[1]
                    self.nbase = x[2]
                elif x[0] == 1:
                    self.x = x[1]
                    self.fx = sum([i * j for i,j in zip(self.x,self.c)])[0]
                    print("Solução ótima:",self.x,"\nValor ótimo da solução:",self.fx)
                elif x[0] == 2:
                    print("Problema não possui solução ótima")