import csv
import numpy as np

np.set_printoptions(precision=4)

matrix_points = []
matrix_polution = []

global points_solution #solução do sistema de pontos amostrais 23,7,30,-50,-12

polution_solution = [0,0,0,0,0] #solução do sistema poluição 22,17,11,72,4

try_pts = [] #chute inicial para os pontos amostrais
try_plt = [] #chute inicial para a poluição

n = 5 #tamanho da matriz

e = 10e-4

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

### calcular chute inicial das matrizes ###
def first_values():
    for i in range(n):
        if matrix_points[i][i] == 0:
            try_pts.append(0)
        else:
            try_pts.append(matrix_points[i][n]/matrix_points[i][i])
        if matrix_polution[i][i] == 0:
            try_plt.append(0)
        else:
            try_plt.append(matrix_polution[i][n]/matrix_polution[i][i])

def sort_lists():

    global points_solution
    points_solution = decLU(matrix_points)[3]
    #print(points_solution)
    for i in range(n):
        for j in range(0, n-i-1):
           if points_solution[j,0] > points_solution[j+1,0]:
               points_solution[j,0], points_solution[j+1,0] = points_solution[j+1,0], points_solution[j,0]

               polution_solution[j], polution_solution[j+1] = polution_solution[j+1], polution_solution[j]

### Função para dados de saída ### 
def saida(): 

    print('-------------------------------------------------')
    print ('|\tPontos  \t|\tPoluição\t|')
    print('-------------------------------------------------')
    for i in range(0,n):

        print('|\t{:.4f}  \t|\t{:.4f}\t\t|'.format(points_solution[i,0], polution_solution[i]))

    print('-------------------------------------------------')
    

###    Decomposição LU    ###
def decLU(A):
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

    x = np.linalg.solve(newU,y)

    return L, newU, y, x


###    Lagrange    ###
def calculates_Lk(x, k):
    lk = 1
    
    for i in range(0,n):
        if i == k:
            continue

        lk = lk * (x - points_solution[i,0])/(points_solution[k,0]-points_solution[i,0])


    return lk

def polinomio(x):
    result = 0

    for i in range (0,n):
        result += polution_solution[i]*calculates_Lk(x,i)

    return result

def main(): 
    
    op = 0
    first_values()
    sort_lists()

    #print(points_solution)

    while True:
        print('------------------------------------------')
        print('1 - Ver tabela (pontos e poluição)')
        print('2 - Resolve o S1 (Decomposição LU)')
        print('3 - Resolve o S2 (Gauss-Seidel)')
        print('4 - Calcula poluição (Interpolação Lagrange)')
        print('5 - Sair')
        op = int(input('Digite a opção: '))

        if op == 1:
            saida()
        elif op == 2:
            ans = decLU(matrix_points)
            print('###    L >> \n{}\n'.format(ans[0]))
            print('###    U >> \n{}\n'.format(ans[1]))
            #print('###    y >> \n{}\n'.format(ans[2]))
            #print('###    x >> \n{}\n'.format(ans[3]))
        elif op == 3:
            print('nothing to do yet')
        elif op == 4:
            x = int(input('Valor do ponto: '))
            print('O grau de poluição é de {:.4f}.\n'.format(polinomio(x)))

        elif op == 5:

            print('Bye bye! :)')
            break
        else: 
            print('Inválido!')


if __name__ == "__main__":
    main()
