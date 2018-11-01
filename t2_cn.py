import csv

matrix = []

#ler matriz
with open('matrix.csv', newline='') as csvfile: 
    line_read = csv.reader(csvfile, delimiter=',')

    for line in line_read:
        matrix.append(line)

