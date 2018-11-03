import csv

matrix_points = []
matrix_polution = []
s1 = [23,7,30,-50,-12]
s2 = [22,17,11,72,4]
n = 5 

#ler matriz
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
        
def sort_lists():
    for i in range(0,n):
        for j in range(0, n-i-1):
           if s1[j] > s1[j+1]:
               s1[j], s1[j+1] = s1[j+1], s1[j]
               s2[j], s2[j+1] = s2[j+1], s2[j]


### Função para dados de saída ### 
def saida(): 

    sort_lists()    

    print('-------------------------------------------------')
    print ('|\tPontos  \t|\tPoluição\t|')
    print('-------------------------------------------------')
    for i in range(0,n):
        print('|\t{:.4f}  \t|\t{:.4f}\t\t|'.format(s1[i], s2[i]))

    print('-------------------------------------------------')


def lineMultiply(l, n):
    new_line = []
    for i in range(len(l)):
        new_line.append(l[i]*n)

    return new_line


def linesSum(la, lb):
    new_line = []
    for i in range(len(la)):
        new_line.append(la[i] + lb[i])
    
    return new_line


def fatLU(A):
    LU = A
    # b
    b = A[n-1]
    # Operações Elementares
    elementarOps = []

    # U
    for i in range(n):
        for j in range(n-1):
            if(i == j and LU[j][i] != 0):
                for k in range(j+1, n):
                    elementarOps.append(-(LU[k][i])/LU[j][i])
                    multipliedLine = lineMultiply(LU[i], elementarOps[-1])
                    LU[k] = linesSum(LU[k], multipliedLine)
    U = LU


    # X
    X = []
    mult = 0
    countX = n
    while(countX != 0):
        U[n][]


    # L
    line = []
    L = []
    count = 0
    for i in range(n):
        L.append(line)
    for j in range(n):
        L[i].append(0)

    



    print('MATRIZ L\n\n')
    for i in range(n):
        print(L[i])
        print()
    print('MATRIZ U\n\n')
    for i in range(n):
        print(U[i])
        print()


    
    #print(elementarOps)


###    Lagrange    ###
def calculates_Lk(x, k):
    lk = 1
    
    for i in range(0,n):
        if i == k:
            continue
        lk = lk * (x - s1[i])/(s1[k]-s1[i])

    return lk

def polinomio(x):
    result = 0

    for i in range (0,n):
        result += s2[i]*calculates_Lk(x,i)

    return result

def main(): 
    
    op = 0

    while True:
        print('------------------------------------------')
        print('1 - Ver tabela (pontos e poluição)')
        print('2 - Resolve A1 (Fatoração LU)')
        print('3 - Resolve A2 (Gauss-Seidel)')
        print('4 - Calcula poluição (Interpolação Lagrange)')
        print('5 - Sair')
        op = int(input('Digite a opção: '))

        if op == 1:
            saida()
        elif op == 2:
            fatLU(matrix_points)
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