import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Определяем исходные данные ---
# Исходная функция
def func(x):
    return np.cos(x) + x

# Узлы для варианта а) (равномерная сетка)
nodes_x_a = np.array([0, np.pi/6, 2*np.pi/6, 3*np.pi/6])
nodes_y_a = func(nodes_x_a)

# Узлы для варианта б) (неравномерная сетка)
nodes_x_b = np.array([0, np.pi/6, np.pi/4, np.pi/2])
nodes_y_b = func(nodes_x_b)

# --- Читаем данные из CSV-файлов ---
try:
    data_a = pd.read_csv('./tasks/3/3.1/plot_data_a.csv')
    data_b = pd.read_csv('./tasks/3/3.1/plot_data_b.csv')
except FileNotFoundError:
    print("Ошибка: CSV-файлы не найдены. Убедитесь, что вы запустили измененный C++ код.")
    exit()

# --- Строим два графика рядом ---
# Создаем фигуру с двумя областями для рисования (1 ряд, 2 колонки)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
fig.suptitle('Сравнение интерполяции для разных наборов узлов', fontsize=16)

# --- График 1: Вариант а) Равномерная сетка ---
ax1.plot(data_a['x'], data_a['actual_f(x)'], label='Исходная функция f(x)', linewidth=2)
ax1.plot(data_a['x'], data_a['polynomial_p(x)'], label='Многочлен P(x)', linestyle='--')
ax1.plot(nodes_x_a, nodes_y_a, 'ro', label='Узлы интерполяции')
ax1.set_title('а) Равномерная сетка')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.legend()

# --- График 2: Вариант б) Неравномерная сетка ---
ax2.plot(data_b['x'], data_b['actual_f(x)'], label='Исходная функция f(x)', linewidth=2)
ax2.plot(data_b['x'], data_b['polynomial_p(x)'], label='Многочлен P(x)', linestyle='--')
ax2.plot(nodes_x_b, nodes_y_b, 'ro', label='Узлы интерполяции')
ax2.set_title('б) Неравномерная сетка')
ax2.set_xlabel('x')
ax2.grid(True)
ax2.legend()

# Показать результат
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()