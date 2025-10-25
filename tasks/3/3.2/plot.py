import numpy as np
import matplotlib.pyplot as plt

# --- ШАГ 1: Вставьте сюда данные из вывода вашей C++ программы ---

# Узлы интерполяции
x_nodes = [0.0, 1.0, 2.0, 3.0, 4.0]
f_nodes = [1.0, 1.5403, 1.5839, 2.01, 3.3464]

# Точка для вычисления
x_star = 1.5

# Коэффициенты сплайнов [a, b, c, d] для каждого интервала
# Скопируйте их из блока "Шаг 2: Коэффициенты кубических сплайнов"
coefficients = [
    # Интервал 1: [0.0, 1.0]
    [1.00000000, 0.68441071, 0.00000000, -0.14411071],
    # Интервал 2: [1.0, 2.0]
    [1.54030000, 0.25207857, -0.43233214, 0.22385357],
    # Интервал 3: [2.0, 3.0]
    [1.58390000, 0.05897500, 0.23922857, 0.12789643],
    # Интервал 4: [3.0, 4.0]
    [2.01000000, 0.92112143, 0.62291786, -0.20763929]
]

# --- Конец блока для вставки данных ---


# Функция для вычисления значения сплайна в точке
def evaluate_spline_segment(x, interval_index):
    """Вычисляет значение одного куска сплайна."""
    a, b, c, d = coefficients[interval_index]
    x_i = x_nodes[interval_index]
    dx = x - x_i
    return a + b * dx + c * dx**2 + d * dx**3

# Настройка и создание графика
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(12, 8))

# 1. Рисуем исходные узлы
ax.plot(x_nodes, f_nodes, 'o', color='red', markersize=8, zorder=10, label='Узлы интерполяции')

# 2. Рисуем каждый кубический сегмент отдельно
for i in range(len(coefficients)):
    # Создаем 100 точек внутри каждого интервала для гладкой кривой
    x_segment = np.linspace(x_nodes[i], x_nodes[i+1], 100)
    y_segment = evaluate_spline_segment(x_segment, i)
    
    # Используем LaTeX для красивого отображения формул в легенде
    label = f'$S_{i+1}(x)$ на [{x_nodes[i]}, {x_nodes[i+1]}]'
    ax.plot(x_segment, y_segment, '--', label=label, linewidth=2)

# 3. Отмечаем точку X* и результат
y_star = evaluate_spline_segment(x_star, 1) # Мы знаем, что 1.5 в интервале 1 (втором)
ax.plot(x_star, y_star, '*', color='green', markersize=15, zorder=11, label=f'S({x_star}) = {y_star:.6f}')
ax.vlines(x_star, 0, y_star, color='green', linestyle=':', linewidth=1.5)


# 4. Оформление
ax.set_title('Визуализация кубического сплайна', fontsize=16)
ax.set_xlabel('x', fontsize=14)
ax.set_ylabel('f(x)', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True)
# Устанавливаем разумные пределы для осей
ax.set_ylim(bottom=0)

# Показать график
plt.show()