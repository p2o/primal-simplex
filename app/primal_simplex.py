import numpy as np

def simplexFii(matrizA,bases,naoBase,vetorC,vetorB):
    B = matrizA[:, bases]
    N = matrizA[:, naoBase]
    cb = vetorC[bases]
    cn = vetorC[naoBase]
    b = vetorB

    # Passo 1: {cálculo da solução básica}
    x = np.linalg.solve(B,b)

    # 2.1) {cálculo do vetor multiplicador simplex}
    lamb = np.linalg.solve(np.transpose(B),cb)

    # 2.2) {custos relativos}
    hat_c = []
    for i in range(N.shape[1]):
        hat_c.append(cn[i,0] - np.dot(np.transpose(lamb),N[:,i])[0])

    # 2.3) {escolha da variável a entrar na base}
    hat_c_min_i = np.argmin(hat_c)
    entra_na_base = naoBase[hat_c_min_i]

    # Passo 3: {teste de otimalidade}
    if hat_c[hat_c_min_i] >= 0 and all(x>=0):
        x_r = [0] * matrizA.shape[1]
        for i,j in enumerate(bases):
            x_r[j] = x[i][0]
        return [1,x_r] #solução ótima

    # Passo 4: {cálculo da direção simplex}
    y = np.linalg.solve(B,N[:,hat_c_min_i])

    # Passo 5: {determinar passo e variável a sair da base}

    # Se y <= 0, então PARE = VERDADE {problema não tem solução ótima finita}
    if max(y) <= 0:
        return [2,x] #sem solução ótima finita

    # Caso contrário, determine a variável a sair da base
    hat_epi = []
    hat_epi_i = []
    for i in range(len(y)):
        if y[i] > 0:
            hat_epi.append(x[i,0] / y[i])
            hat_epi_i.append(i)

    hat_epi_min_i = hat_epi_i[np.argmin(hat_epi)]
    sai_da_base = bases[hat_epi_min_i]

    # Passo 6: {atualização: nova partição básica}
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
        A_columns = [i for i in range(self.A.shape[1])]
        self.base = [-1] * tamanho_base
        self.nbase = []
        id_columns_found = []
        matriz_B_to_find = np.identity(tamanho_base)

        # Verifica se a matriz A possui base Identidade para usar como partição básica factível inicial
        for col_B in range(tamanho_base):
            for col_A in A_columns:
                if all(self.A[:,col_A] == matriz_B_to_find[:,col_B]):
                    self.base[len(id_columns_found)] = col_A
                    id_columns_found.append(col_B)
                    A_columns.remove(col_A)
                    break
            if len(id_columns_found) >= tamanho_base:
                self.nbase = [i for i in range(self.A.shape[1]) if i not in self.base]
                break

        # Se a quantidade de colunas encontradas para base é menor que o tamanho necessário para a base então:
        # é criado um problema auxiliar para encontrar uma base inicial factível utilizando a fase II
        if len(id_columns_found) < tamanho_base:
            columns_to_add = np.asarray([matriz_B_to_find[i] for i in range(tamanho_base) if i not in id_columns_found]).T
            A_aux = np.concatenate((self.A,columns_to_add), axis=1)
            c_aux = np.asarray([0] * self.A.shape[1] + [1] * columns_to_add.shape[1]).reshape(A_aux.shape[1],1)
            base_i = len(id_columns_found)
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
        self.solucao = [0,""]
        if not hasattr(self, "base"):
            self.solveFi()
        if self.A.shape[1] < len(self.base) + len(self.nbase):
            self.solucao = [3,"Problema infactivel"]
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
                    self.solucao = [1,"Solução ótima: " + str(self.x) + " Valor ótimo da solução: " + str(self.fx)]
                elif x[0] == 2:
                    self.solucao = [2,"Problema não possui solução ótima"]

        print(self.solucao[1])

        return self.solucao
