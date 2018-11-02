import csv

matrix_points = []
matrix_polution = []
points_solution = [23,7,30,-50,-12] #solução do sistema de pontos amostrais 23,7,30,-50,-12
polution_solution = [0,0,0,0,0] #solução do sistema poluição 22,17,11,72,4

try_pts = [] #chute inicial para os pontos amostrais
try_plt = [] #chute inicial para a poluição
n = 5 
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
    for i in range(n):
        for j in range(0, n-i-1):
           if points_solution[j] > points_solution[j+1]:
               points_solution[j], points_solution[j+1] = points_solution[j+1], points_solution[j]
               polution_solution[j], polution_solution[j+1] = polution_solution[j+1], polution_solution[j]

### Função para dados de saída ### 
def saida(): 

    gauss_Seidel()

    print('-------------------------------------------------')
    print ('|\tPontos  \t|\tPoluição\t|')
    print('-------------------------------------------------')
    for i in range(0,n):
        print('|\t{:.4f}  \t|\t{:.4f}\t\t|'.format(points_solution[i], polution_solution[i]))

    print('-------------------------------------------------')

### Erro ###
def calculates_Error(x,m): # Recebe o vetor anterior e a indicação da matriz usada (1 para A1 e 2 para A2)
    maxDiff = 0
    diff = 0
    maxSol = 0
    for i in range(n):
        if m == 1:
            diff = abs(try_pts[i] - x[i])
        elif m == 2:
            diff = abs(try_plt[i] - x[i])
            print('diff {} {} {}'.format(x[i], try_plt[i], diff))
        if diff > maxDiff:
            maxDiff = diff
        if x[i] > maxSol:
            maxSol = x[i]
    if maxSol == 0:
        return maxDiff   
    return maxDiff/maxSol


### Método Gauss-Seidel ###
def gauss_Seidel():
    k = 0 
    aux = 0
    xAnt = [0 for i in range(n)]
    while(calculates_Error(xAnt, 2)):
        for i in range(0,n): 
            print(matrix_polution[i][n-1])
            aux = matrix_polution[i][n-1]
            for j in range(n):
                if i == j:
                    continue
                aux = aux - matrix_polution[i][j]*try_plt[j]
            xAnt[i] = try_plt[i]
            if matrix_polution[i][i] == 0:
                try_plt[i] = aux
            else:
                try_plt[i] = aux/matrix_polution[i][i]



###    Lagrange    ###
def calculates_Lk(x, k):
    lk = 1
    
    for i in range(0,n):
        if i == k:
            continue
        lk = lk * (x - points_solution[i])/(points_solution[k]-points_solution[i])

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

    while True:
        print('------------------------------------------')
        print('1 - Ver tabela (pontos e poluição)')
        print('2 - Fatoração LU (Processo)')
        print('3 - Gauss-Seidel (Processo)')
        print('4 - Calcula poluição (Interpolação Lagrange)')
        print('5 - Sair')
        op = int(input('Digite a opção: '))

        if op == 1:
            saida()
        elif op == 2:
            x = int(input('Valor do ponto: '))
            print('O grau de poluição é de {:.4f}.\n'.format(polinomio(x)))
        elif op == 3:
            gauss_Seidel()
        elif op == 5: 
            print('Bye bye! :)')
            break
        else: 
            print('Inválido!')


if __name__ == "__main__":
    main()