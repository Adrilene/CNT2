import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange # apenas para printar o polinômio de forma visível e plotar o gráfico


n = 5 #tamanho da matriz
np.set_printoptions(precision=4)
matrix_points = [] # matriz A
matrix_polution = [] # matriz B 
points_solution = np.zeros([n],dtype=float)  # Armazena solução do S1
#print(points_solution)
polution_solution = np.zeros([n],dtype=float) # Armazena a solução do S2
#print(polution_solution)


### ler matrizes ###
with open('matrix.csv', newline='') as csvfile: 
    line_read = csv.reader(csvfile, delimiter=',')
    m = 1
    for line in line_read:
        if line == []:
            m = 2
            continue
        if m == 1:
            matrix_points.append(list(map(float,line)))
        elif m == 2:
            matrix_polution.append(list(map(float,line)))


###   Verificação   ###
def verificacao():
    #Verifica se os sistemas já foram solucionados,
    #caso contrário encerra a função informando qual não foi resolvido.
    if(max(points_solution) == 0 and max(polution_solution) == 0):
        print('####    ATENÇÃO    ####\nResolva o S1 e o S2 primeiro!\n')
        return 0
    elif(max(points_solution) == 0):
        print('####    ATENÇÃO    ####\nResolva o S1 primeiro!\n')
        return 0
    elif(max(polution_solution) == 0):
        print('####    ATENÇÃO    ####\nResolva o S2 primeiro!\n')
        return 0
    else:
        return 1



###  Ordenar Listas ###
def sort_lists():
    
    if not verificacao():
        return 0
    
    #ordena a points_solution
    #o mesmo swap feito nela é aplicado a polution_solution, para que os resultados sejam respectivos
    for i in range(n):
        for j in range(n-i-1):
            if points_solution[j] > points_solution[j+1]:
                points_solution[j], points_solution[j+1] = points_solution[j+1], points_solution[j]
                polution_solution[j], polution_solution[j+1] = polution_solution[j+1], polution_solution[j]
        
    return 1

### Função para dados de saída ### 
def saida(): 

    if(sort_lists()):
        print('-------------------------------------------------')
        print ('|\tPontos  \t|\tPoluição\t|')
        for i in range(0,n):

            print('|\t{:.4f}  \t|\t{:.4f}\t\t|'.format(points_solution[i], polution_solution[i]))


        print('-------------------------------------------------')
    else:
        pass
    

###    Decomposição LU    ###
def decLU():

    global points_solution

    U = np.array(matrix_points)
    L = np.eye(n)
    
    for j in range(n-1):
        for i in range(j+1, n):
            L[i,j] = U[i,j]/U[j,j]
            for k in range(j+1, n):
                U[i,k] = U[i,k] - (L[i,j] * U[j,k])
            U[i,j] = 0

    b = np.delete(U, np.s_[0:n], 1)
    
    newU = np.delete(U, np.s_[n], 1)

    y = np.linalg.solve(L,b)
    yn = np.reshape(y,n) # array

    x = np.linalg.solve(newU,y)
    xn = np.reshape(x,n) # array

    points_solution = xn

    return L, newU, yn, xn


### Método Gauss-Seidel ###

def sassenfeld():
    b = [1 for i in range(n)]
    b[0] = 0

    for i in range(1,n): 
        b[0] = b[0] + matrix_polution[0][i]

    b[0] = b[0]/matrix_polution[0][0]

    for i in range(1,n):
        for j in range(n):
            if j == i:
                continue
            for k in range(j-1,-1,-1):
                b[i] = matrix_polution[i][j] * b[k]
            b[i] = b[i]/matrix_polution[j][j]
    
    if max(b) >= 1:
        return False
    return True

### Calcula Chute inicial para Matriz de Poluição (B) ###

def chuteInicial():
    x = []
    for i in range(n): 
        if matrix_polution[i][0] != 0:
            x.append(matrix_polution[i][n]/matrix_polution[i][0])
        else: 
            x.append(0)

    return x

def calculaErro(previous):
    
    maxDiff = 0 #Armazenará o máximo da diferença entre os resultados
    diff = 0 #Auxiliar
    maxSol = max(polution_solution) #Armazena o máximo da solução atual

    if maxSol == 0: #Se o máximo da solução for zero, retorna 1, pois não pode dividir por zero.
        return 1
    
    for i in range(n):
        diff = abs(polution_solution[i] - previous[i])
        
        if diff > maxDiff:
            maxDiff = diff
    
    return maxDiff/maxSol


#Recebe o valor de erro, número máximo de iterações e uma flag para indicar se irá mostrar as iterações
def gaussSeidel(erro, maxInt, show): 
    
    previousX = chuteInicial() #Recebe o chute inicial
    global polution_solution

    k = 0 
    
    while(k<maxInt and calculaErro(previousX)>erro):
        for i in range(n):
            previousX[i] = polution_solution[i]
            soma = 0 

            for j in range(n):
                if j == i: 
                    continue
                soma = soma + matrix_polution[i][j] * polution_solution[j]

            polution_solution[i] = (matrix_polution[i][n]-soma)/matrix_polution[i][i]

        if show == 1:
            print('Solução na iteração {}: {}'.format(k,polution_solution))
            print('Erro: {}'.format(calculaErro(previousX)))
            print('\n')

        k = k+1

###    Lagrange    ###

def calculates_Lk(x, k):
    lk = 1
    
    for i in range(0,n):
        if i == k:
            continue

        lk = lk * (x - points_solution[i])/(points_solution[k]-points_solution[i])

    return lk

def result_Polinomio(x): 
    result = 0

    for i in range (n):
        result = result + polution_solution[i]*calculates_Lk(x,i)
    
    return result


def plot_results(lx,ly):
    sort_lists()
    
    x = points_solution
    y = polution_solution
    city = (lx,ly)

    plt.figure()
    plt.title('Pontos Amostrais  X  Poluição\nInterpolação Linear(Lagrange)')

    plt.xlabel('Pontos Amostrais(PA)')
    plt.ylabel('Poluição(P)')
    plt.plot(x, y, 'ro', label='PA-P')
    plt.plot(x, y, 'k--', label='Interpolação')
    plt.plot(city[0],city[1], 'bo', label='Cidade')

    plt.legend()
    plt.show()
    

def main(): 
    op = str(0)

    #print(points_solution)

    while True:

        if op == str(0):
            print('------------------------------------------')
            print('1 - Resolve o S1 (Decomposição LU)')
            print('2 - Resolve o S2 (Gauss-Seidel)')
            print('3 - Ver tabela (pontos e poluição)')
            print('4 - Calcula poluição (Interpolação Lagrange)')
            print('5 - Sair')
            op = str(input('Digite a opção: '))
            print('\n')
        elif op == str(1):
            op = str(0)
            ans = decLU()
            print('###    L >> \n{}\n'.format(ans[0]))
            print('###    U >> \n{}\n'.format(ans[1]))
            print('###    y >> \n{}\n'.format(ans[2]))
            print('###    x >> \n{}\n'.format(ans[3]))
        elif op == str(2):
            op = str(0)
            #testa se critério de convergência é satisfeito
            if sassenfeld():
                try:
                    e = float(input('Erro: '))
                except ValueError:
                    print('####    ATENÇÃO    ####\ntype(Erro) = float\n')
                    continue
                try:
                    maxInt = int(input('Número máximo de Interações: '))
                except ValueError:
                    print('####    ATENÇÃO    ####\ntype(Num) = int\n')
                    continue

                gaussSeidel(e,maxInt,1)
            else:
                print('A matriz não irá convergir.\n')
        elif op == str(3):
            op = str(0)
            saida()
        elif op == str(4):
            op = str(0)
            if not verificacao():
                continue

            print('p(x) = \n{}'.format(lagrange(points_solution,polution_solution)))
            x = float(input('Valor do ponto(x): '))
            y = result_Polinomio(x)
            print('O grau de poluição é de {:.4f}.\n'.format(y))
            opt = 0
            while True:
                if opt == 's':
                    plot_results(x,y)
                    print()
                    break
                elif opt == 'n':
                    print()
                    break
                else:
                    opt = str(input('Deseja plotar o gráfico?(s/n): '))
        elif op == str(5):
            print('Bye bye! :)')
            break
        else: 
            op = str(0)
            print('####    ATENÇÃO    ####\nOpção Inválida!\n')


if __name__ == "__main__":
    main()
