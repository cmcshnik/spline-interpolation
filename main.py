'''
Кольцов К.Е., ВМК МГУ, группа 209
Интерполирование функций 

Для отрисовки графиков используется библиотека matplotlib.
Все ключевые параметры программы настраиваются через глобальные 
переменные, определенные ниже.
'''

import matplotlib.pyplot as plt
import math

# Ключевые параметры 
FUNC_NUM = 1           # номер интерполируемой функции (1 или 2)
INTERPOLATOR_NUM = 1   # номер интерполятора 
                       # (1 - лагранж, 2 - сплайны)

LEFT_BOUNDARY = 0  # Левая граница отрезка
RIGHT_BOUNDARY = 2  # Правая граница отрезка
INTERPOLATION_NODES_COUNT = 17  # Количество узлов интерполирования
                                # для отрисовки на графике


CONVERGENCE_DO_CHECK = 0  # Делать ли проверку на сходимость (1 - да, 0 - нет)
CONVERGENCE_NODES_LIMIT = 80  # Макс. число узлов интерполяции
CONVERGENCE_POINTS_AMOUNT = 1001  # Количество точек в сетке 
                                  # для сравнения значений функций
CONVERGENCE_GRAPH_FILE_NAME = "conv.png"  # Название файла, в который будет записан
                                          # график изменения разницы значений функций в точках
CONVERGENCE_DO_PRINT = 1  # Печатать ли в консоль информацию о сходимости (1 - да, 0 - нет)
                                

GRAPH_FILE_NAME = "graph.png"  # Название файла, в который будут записаны графики функций
EPSILON = 1e-10  # Величина, ниже которой разность между двумя числами 
                 # считается малой и числа признаются равными
STEP_FOR_PLOTTING = 0.001  # То, с каким шагом будут 
                           # высчитываться значения функций 
                           # для построения графика 

# Первая функция в упрощенном виде
def Func1(x):
    return 63/8 * (x**5) - 35/4 * (x**3) + 15/8 * (x)

def Func2(x):
    return abs(math.sin(4 * x)) * math.exp(2 * x)

# Функция, считающая и возвращающая список узлов равномерной сетки
def CreateNodes(amount):
    x_i = LEFT_BOUNDARY
    res = [x_i]
    step = (RIGHT_BOUNDARY - LEFT_BOUNDARY) / amount

    while (x_i <= RIGHT_BOUNDARY):
        x_i += step
        res.append(x_i)

    return res

# Отрисовка функции с заданным шагом
def PlottingFunction(func, func_name, color):
    arr = CreateNodes((RIGHT_BOUNDARY - LEFT_BOUNDARY) / STEP_FOR_PLOTTING)
    y = [func(i) for i in arr] 

    plt.plot(arr, y, label=func_name, color=color)  
    plt.legend() 

# Функция, считающая значения функции func в узлах nodes
def CreateValues(func, nodes):
    res = []

    for node in nodes:
        res.append(func(node))

    return res


def RunThroughMethod(a, b, c, d):
    '''
    Решение трёхдиагональной системы методом прогонки.
    
    a: поддиагональ (длина N-1)
    b: главная диагональ (длина N)
    c: наддиагональ (длина N-1)
    d: правая часть (длина N)
    
    Возвращает: x — решение системы (длина N)
    '''
    
    n = len(b)
    
    c_prime = [0] * (n - 1)  
    d_prime = [0] * n        
    
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n):
        denom = b[i] - a[i - 1] * c_prime[i - 1]
        if i < n - 1:
            c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / denom
    
    x = [0] * n
    x[-1] = d_prime[-1]
    
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    
    return x

class SplineInterpolant:
    '''
    Интерполяция кубическимим сплайнами по точкам и значениям.
    При создании объекта происходит вычисление коэффициентов полиномов.
    '''

    def __init__(self, nodes, values):
        self.nodes = nodes
        self.values = values
        self.N = len(nodes) - 1
        nodes_amount = len(nodes)

        # ГРАНИЧНЫЕ ЗНАЧЕНИЯ КОЭФФИЦИЕНТА с НУЛЕВЫЕ
        c0 = 0
        cN = 0
        
        '''
            Все массивы имеют избыточный размер для реализации 
            индексации ровно как в формулах. 
            Во всех циклах и срезах, где фигурирует N (self.N),
            второй предел увеличен на единицу
            в связи с оссобенностями работы python
        '''
        self.a = [0 for _ in range(nodes_amount)]
        self.b = [0 for _ in range(nodes_amount)]
        self.c = []
        self.d = [0 for _ in range(nodes_amount)]
        h = [0 for _ in range(nodes_amount)]

        # Подсчёт коэффициентов а
        for i in range(len(values)):
            self.a[i] = values[i]

        # Подсчёт коэффициентов c
        for i in range(1, self.N + 1):
            h[i] = nodes[i] - nodes[i - 1]
        
        temp_h = [0 for _ in range(nodes_amount)]
        temp_f = [0 for _ in range(nodes_amount)]
        for i in range(1, (self.N - 1) + 1):
            temp_h[i] = 2 * (h[i] + h[i + 1])
            temp_f[i] = (6 * ((values[i + 1] - values[i]) / h[i + 1] - (values[i] - values[i - 1]) / h[i]))

        temp_f[1] -= h[1] * c0
        temp_f[self.N-1] -= h[self.N] * cN 
        
        self.c = RunThroughMethod(h[2:(self.N-1)+1], temp_h[1:(self.N-1)+1], h[2:(self.N - 1) + 1], temp_f[1:(self.N-1)+1])
        self.c.insert(0, c0)
        self.c.append(cN)

        # Подсчёт коэффициентов d
        for i in range(1, self.N + 1):
            self.d[i] = (self.c[i] - self.c[i - 1]) / h[i]

        # Подсчёт коэффициентов b
        for i in range(1, self.N + 1):
            self.b[i] = h[i] * self.c[i] / 2 - h[i]*h[i] * self.d[i] / 6 + (values[i] - values[i - 1]) / h[i]
        
    # Подсчёт значения для случайного х
    def Calc(self, x):
        for i in range(1, self.N + 1):
            if (self.nodes[i - 1] <= x <= self.nodes[i]):
                dx = x - self.nodes[i]
                res = self.a[i] + self.b[i] * dx + self.c[i] * dx*dx / 2 + self.d[i] * dx*dx*dx / 6
                return res


class LagrangeInterpolant:
    def __init__(self, nodes, values):
        self.nodes = nodes
        self.values = values
        self.n = len(nodes)
        
        # Предварительный расчет барицентрических весов
        self.c = [1] * self.n
        for j in range(self.n):
            for k in range(self.n):
                if k != j:
                    self.c[j] /= (self.nodes[j] - self.nodes[k])

    # Подсчёт значения для случайного х
    def Calc(self, x):
        numerator = 0
        denominator = 0

        # Вычисляем числитель и знаменатель барицентрической формулы
        for j in range(self.n):
            # Проверка на совпадение с узлом
            if abs(x - self.nodes[j]) < EPSILON:
                return self.values[j]
    
            term = self.c[j] / (x - self.nodes[j])
            numerator += term * self.values[j]
            denominator += term

        return numerator / denominator


def CheckConvergence(func):
    '''
    Проверка на сходимость
    func          - интерполируемая функция
    num_of_points - то, сколько узлов будет в сетке, на которой
                    происходит сравнение значний функций
    nodes_limit   - лимит узлов интерполирования, >= 2
    '''

    max_deviation = 0
    values = []
    nodes = []
    max_deviations = []
    num_of_points = CONVERGENCE_POINTS_AMOUNT
    nodes_limit = CONVERGENCE_NODES_LIMIT

    def Compare(func1, func2, num_of_points):
        max_deviation = 0
        step = (RIGHT_BOUNDARY - LEFT_BOUNDARY) / num_of_points

        disp = LEFT_BOUNDARY
        while (disp <= RIGHT_BOUNDARY):
            diff = abs(func1(disp) - func2(disp))

            if (diff > max_deviation):
                max_deviation = diff

            disp += step

        return max_deviation

    for node_num in range(3, nodes_limit+1):
        nodes = CreateNodes(node_num)
        values = CreateValues(func, nodes)
        
        if (INTERPOLATOR_NUM == 1):
            P = LagrangeInterpolant(nodes, values)
        elif (INTERPOLATOR_NUM == 2):
            P = SplineInterpolant(nodes, values)

        max_deviation = Compare(func, P.Calc, num_of_points)
        max_deviations.append(max_deviation)

        if (CONVERGENCE_DO_PRINT == 1):
            print(f"Кол-во узлов: {node_num}; Макс. отклонение: {max_deviation}")
    
    plt.figure()
    plt.title(f"График отклонения для функции f{FUNC_NUM}")
    plt.xlabel("Количество узлов")  
    plt.ylabel("Максимальная разница значений") 
    plt.grid(True)  
    plt.plot([x for x in range(3, nodes_limit+1)], max_deviations) 
    plt.savefig(CONVERGENCE_GRAPH_FILE_NAME) 

def main(): 
    plt.figure()
    plt.xlabel("x")  
    plt.ylabel("y") 
    plt.title(f"Число узлов интерполяции: {INTERPOLATION_NODES_COUNT}")
    plt.grid(True)  

    if (FUNC_NUM == 1):
        nodes = CreateNodes(INTERPOLATION_NODES_COUNT)
        values = CreateValues(Func1, nodes)

        PlottingFunction(Func1, "Function 1", "blue")

        if (INTERPOLATOR_NUM == 1):
            P = LagrangeInterpolant(nodes, values)
            PlottingFunction(P.Calc, "Lagrange", "red")
        elif (INTERPOLATOR_NUM == 2):
            P = SplineInterpolant(nodes, values)
            PlottingFunction(P.Calc, "Spline", "red")

        plt.savefig(GRAPH_FILE_NAME)

        if (CONVERGENCE_DO_CHECK):
            CheckConvergence(Func1)
    elif (FUNC_NUM == 2):
        nodes = CreateNodes(INTERPOLATION_NODES_COUNT)
        values = CreateValues(Func2, nodes)

        PlottingFunction(Func2, "Function 2", "blue")

        if (INTERPOLATOR_NUM == 1):
            P = LagrangeInterpolant(nodes, values)
            PlottingFunction(P.Calc, "Lagrange", "red")
        elif (INTERPOLATOR_NUM == 2):
            P = SplineInterpolant(nodes, values)
            PlottingFunction(P.Calc, "Spline", "red")

        plt.savefig(GRAPH_FILE_NAME)

        if (CONVERGENCE_DO_CHECK):
            CheckConvergence(Func2)

if __name__ == "__main__":
    main()
