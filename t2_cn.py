import csv

n = 5 

matrix_points = []
matrix_polution = []
points_solution = [23,7,30,-50,-12] #solução do sistema de pontos amostrais 23,7,30,-50,-12
polution_solution = [22,17,11,72,4] #solução do sistema poluição 22,17,11,72,4

try_pts = [] #chute inicial para os pontos amostrais
try_plt = [] #chute inicial para a poluição

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

###  Ordenar Listas ###
def sort_lists():
    for i in range(n):
        for j in range(0, n-i-1):
           if points_solution[j] > points_solution[j+1]:
               points_solution[j], points_solution[j+1] = points_solution[j+1], points_solution[j]
               polution_solution[j], polution_solution[j+1] = polution_solution[j+1], polution_solution[j]

### Função para dados de saída ### 
def saida(): 

    sort_lists()

    print('-------------------------------------------------')
    print ('|\tPontos  \t|\tPoluição\t|')
    for i in range(0,n):
        print('-------------------------------------------------')
        print('|\t{:.4f}  \t|\t{:.4f}\t\t|'.format(points_solution[i], polution_solution[i]))

    print('-------------------------------------------------')

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

    for i in range (0,n):
        result += polution_solution[i]*calculates_Lk(x,i)

    return result

### Fim Lagrange ### 

def main(): 
    
    op = 0

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
            pass
        elif op == 3:
            gauss_Seidel()
        elif op == 4:
            x = int(input("Digite o ponto: "))
            print("A poluição é de: {:.4f}.".format(result_Polinomio(x)))
        elif op == 5: 
            print('Bye bye! :)')
            break
        else: 
            print('Inválido!')


if __name__ == "__main__":
    main()