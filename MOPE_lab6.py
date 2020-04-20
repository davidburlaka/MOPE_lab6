from numpy.linalg import solve
from _pydecimal import Decimal
from scipy.stats import f
from scipy.stats import t
from random import randrange
from math import sqrt
from math import fabs as fab

dep, undep = 0, 0
m, d = 3, 0
p = 0.95
N = 15

x1_min, x1_max = 15, 45
x2_min, x2_max = 30, 80
x3_min, x3_max = 15, 45
x01 = (x1_max + x1_min) / 2
x02 = (x2_max + x2_min) / 2
x03 = (x3_max + x3_min) / 2
delta_x1 = x1_max - x01
delta_x2 = x2_max - x02
delta_x3 = x3_max - x03

matrix_pfe = [
    [-1, -1, -1, +1, +1, +1, -1, +1, +1, +1],
    [-1, -1, +1, +1, -1, -1, +1, +1, +1, +1],
    [-1, +1, -1, -1, +1, -1, +1, +1, +1, +1],
    [-1, +1, +1, -1, -1, +1, -1, +1, +1, +1],
    [+1, -1, -1, -1, -1, +1, +1, +1, +1, +1],
    [+1, -1, +1, -1, +1, -1, -1, +1, +1, +1],
    [+1, +1, -1, +1, -1, -1, -1, +1, +1, +1],
    [+1, +1, +1, +1, +1, +1, +1, +1, +1, +1],
    [-1.73, 0, 0, 0, 0, 0, 0, 2.9929, 0, 0],
    [+1.73, 0, 0, 0, 0, 0, 0, 2.9929, 0, 0],
    [0, -1.73, 0, 0, 0, 0, 0, 0, 2.9929, 0],
    [0, +1.73, 0, 0, 0, 0, 0, 0, 2.9929, 0],
    [0, 0, -1.73, 0, 0, 0, 0, 0, 0, 2.9929],
    [0, 0, +1.73, 0, 0, 0, 0, 0, 0, 2.9929],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

class Critical_values:
    @staticmethod
    def get_cohren_value(size_of_selections, qty_of_selections, significance):
        size_of_selections += 1
        partResult1 = significance / (size_of_selections - 1)
        params = [partResult1, qty_of_selections, (size_of_selections - 1 - 1) * qty_of_selections]
        fisher = f.isf(*params)
        result = fisher / (fisher + (size_of_selections - 1 - 1))
        return Decimal(result).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def get_student_value(f3, significance):
        return Decimal(abs(t.ppf(significance / 2, f3))).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def get_fisher_value(f3, f4, significance):
        return Decimal(abs(f.isf(significance, f4, f3))).quantize(Decimal('.0001')).__float__()


def generate_matrix():
    def f(X1, X2, X3):
        y = 0.3 + 4.1 * X1 + 2.8 * X2 + 7.8 * X3 + 1.4 * X1 * X1 + 0.2 * X2 * X2 + 2.4 * X3 * X3 + 9.7 * X1 * X2 + \
            0.6 * X1 * X3 + 4.4 * X2 * X3 + 3.4 * X1 * X2 * X3 + randrange(0, 10) - 5
        return y

    matrix_y = [[f(matrix_x[j][0], matrix_x[j][1], matrix_x[j][2]) for i in range(m)] for j in range(N)]
    return matrix_y


def x(l1, l2, l3):
    x_1 = l1 * delta_x1 + x01
    x_2 = l2 * delta_x2 + x02
    x_3 = l3 * delta_x3 + x03
    return [x_1, x_2, x_3]


def find_average(lst, orientation):
    average = []
    if orientation == 1:  # Середнє значення по рядку
        for rows in range(len(lst)):
            average.append(sum(lst[rows]) / len(lst[rows]))
    else:  # Середнє значення по колонкі
        for column in range(len(lst[0])):
            number_lst = []
            for rows in range(len(lst)):
                number_lst.append(lst[rows][column])
            average.append(sum(number_lst) / len(number_lst))
    return average


def a(first, second):  # first = 1, second = 2 : пошук а12
    need_a = 0
    for j in range(N):
        need_a += matrix_x[j][first - 1] * matrix_x[j][second - 1] / N
    return need_a


def find_known_a(number):
    a = 0
    for j in range(N):
        a += average_y[j] * matrix_x[j][number - 1] / 15
    return a


def check_result(b_lst, k):
    y_i = b_lst[0] + b_lst[1] * matrix[k][0] + b_lst[2] * matrix[k][1] + b_lst[3] * matrix[k][2] + \
          b_lst[4] * matrix[k][3] + b_lst[5] * matrix[k][4] + b_lst[6] * matrix[k][5] + b_lst[7] * matrix[k][6] + \
          b_lst[8] * matrix[k][7] + b_lst[9] * matrix[k][8] + b_lst[10] * matrix[k][9]
    return y_i


def student(b_lst, number_x=10):
    disp_b = sqrt(disp_b2)
    for column in range(number_x + 1):
        tp = 0
        tt = Critical_values.get_student_value(f3, q)
        for row in range(N):
            if column == 0:
                tp += average_y[row] / N
            else:
                tp += average_y[row] * matrix_pfe[row][column - 1]
        if fab(tp / disp_b) < tt:
            b_lst[column] = 0
    return b_lst


def fisher():
    disp_ad = 0
    f4 = N - d
    for row in range(len(average_y)):
        disp_ad += (m * (average_y[row] - check_result(student_lst, row))) / (N - d)
    Fp = disp_ad / disp_b2
    Ft = Critical_values.get_fisher_value(f3, f4, q)
    return Fp < Ft

for i in range(100):
    matrix_x = [[] for x in range(N)]
    for i in range(len(matrix_x)):
        if i < 8:
            x_1 = x1_min if matrix_pfe[i][0] == -1 else x1_max
            x_2 = x2_min if matrix_pfe[i][1] == -1 else x2_max
            x_3 = x3_min if matrix_pfe[i][2] == -1 else x3_max
        else:
            x_lst = x(matrix_pfe[i][0], matrix_pfe[i][1], matrix_pfe[i][2])
            x_1, x_2, x_3 = x_lst
        matrix_x[i] = [x_1, x_2, x_3, x_1 * x_2, x_1 * x_3, x_2 * x_3, x_1 * x_2 * x_3, x_1 ** 2, x_2 ** 2, x_3 ** 2]

    adequacy, homogeneity = False, False
    while not adequacy:
        matrix_y = generate_matrix()
        average_x = find_average(matrix_x, 0)
        average_y = find_average(matrix_y, 1)
        matrix = [(matrix_x[i] + matrix_y[i]) for i in range(N)]
        mx_i = average_x
        my = sum(average_y) / 15

        unknown = [
            [1, mx_i[0], mx_i[1], mx_i[2], mx_i[3], mx_i[4], mx_i[5], mx_i[6], mx_i[7], mx_i[8], mx_i[9]],
            [mx_i[0], a(1, 1), a(1, 2), a(1, 3), a(1, 4), a(1, 5), a(1, 6), a(1, 7), a(1, 8), a(1, 9), a(1, 10)],
            [mx_i[1], a(2, 1), a(2, 2), a(2, 3), a(2, 4), a(2, 5), a(2, 6), a(2, 7), a(2, 8), a(2, 9), a(2, 10)],
            [mx_i[2], a(3, 1), a(3, 2), a(3, 3), a(3, 4), a(3, 5), a(3, 6), a(3, 7), a(3, 8), a(3, 9), a(3, 10)],
            [mx_i[3], a(4, 1), a(4, 2), a(4, 3), a(4, 4), a(4, 5), a(4, 6), a(4, 7), a(4, 8), a(4, 9), a(4, 10)],
            [mx_i[4], a(5, 1), a(5, 2), a(5, 3), a(5, 4), a(5, 5), a(5, 6), a(5, 7), a(5, 8), a(5, 9), a(5, 10)],
            [mx_i[5], a(6, 1), a(6, 2), a(6, 3), a(6, 4), a(6, 5), a(6, 6), a(6, 7), a(6, 8), a(6, 9), a(6, 10)],
            [mx_i[6], a(7, 1), a(7, 2), a(7, 3), a(7, 4), a(7, 5), a(7, 6), a(7, 7), a(7, 8), a(7, 9), a(7, 10)],
            [mx_i[7], a(8, 1), a(8, 2), a(8, 3), a(8, 4), a(8, 5), a(8, 6), a(8, 7), a(8, 8), a(8, 9), a(8, 10)],
            [mx_i[8], a(9, 1), a(9, 2), a(9, 3), a(9, 4), a(9, 5), a(9, 6), a(9, 7), a(9, 8), a(9, 9), a(9, 10)],
            [mx_i[9], a(10, 1), a(10, 2), a(10, 3), a(10, 4), a(10, 5), a(10, 6), a(10, 7), a(10, 8), a(10, 9), a(10, 10)]
        ]
        known = [my, find_known_a(1), find_known_a(2), find_known_a(3), find_known_a(4), find_known_a(5), find_known_a(6),
                find_known_a(7), find_known_a(8), find_known_a(9), find_known_a(10)]

        beta = solve(unknown, known)
        print("Отримане рівняння регресії")
        print("ŷ = {:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
          "+ {:.3f} * Х1Х2X3 + {:.3f} * X11² + {:.3f} * X22² + {:.3f} * X33²\nПеревірка"
          .format(beta[0], beta[1], beta[2], beta[3], beta[4], beta[5], beta[6], beta[7], beta[8], beta[9], beta[10]))
        for i in range(N):
            print("ŷ{} = {:.3f} ≈ {:.3f}".format((i + 1), check_result(beta, i), average_y[i]))

        while not homogeneity:
            print("\n{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}{:^13}"
              .format("X1", "X2", "X3", "X1X2", "X1X3", "X2X3", "X1X2X3", "X1X1", "X2X2", "X3X3", "Y1", "Y2", "Y3"))
            for row in range(N):
                for column in range(len(matrix[0])):
                    print("{:^12.3f}".format(matrix[row][column]), end=' ')
                print()
            disp_y = [0.0 for x in range(N)]
            for i in range(N):
                disp_i = 0
                for j in range(m):
                    disp_i += (matrix_y[i][j] - average_y[i]) ** 2
                disp_y.append(disp_i / (m - 1))
            f1 = m - 1
            f2 = N
            f3 = f1 * f2
            q = 1 - p
            Gp = max(disp_y) / sum(disp_y)
            print("\nПеревірка за Кохреном")
            Gt = Critical_values.get_cohren_value(f2, f1, q)
            if Gt > Gp:
                print("Дисперсія однорідна при рівні значимості {:.2f}\nЗбільшувати m не потрібно.".format(q))
                homogeneity = True
            else:
                print("Дисперсія не однорідна при рівні значимості {:.2f}".format(q))
                m += 1

        disp_b2 = sum(disp_y) / (N * N * m)
        student_lst = list(student(beta))
        print("\nОтримане рівняння регресії з урахуванням критерія Стьюдента")
        print("ŷ = {:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
          "+ {:.3f} * Х1Х2X3 + {:.3f} * X11² + {:.3f} * X22² + {:.3f} * X33²\nПеревірка"
          .format(student_lst[0], student_lst[1], student_lst[2], student_lst[3], student_lst[4], student_lst[5],
                  student_lst[6], student_lst[7], student_lst[8], student_lst[9], student_lst[10]))
        for i in range(N):
            print("ŷ{} = {:.3f} ≈ {:.3f}".format((i + 1), check_result(student_lst, i), average_y[i]))

        print("\nПеревірка за Фішером")
        d = 11 - student_lst.count(0)
        dep += N - d
        undep += d
        if fisher():
            print("Рівняння регресії адекватне оригіналу")
            adequacy = True
        else:
            print("Рівняння регресії неадекватне оригіналу")
print("Кільість значимих: ", dep)
print("Кількість не значимих: ", undep)
