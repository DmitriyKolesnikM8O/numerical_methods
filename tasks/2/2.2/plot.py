import numpy as np
import matplotlib.pyplot as plt

# --- Параметры ---
# Начальное приближение, которое мы определили, посмотрев на график
X1_INIT = 0.42
X2_INIT = 0.97

# --- Определяем наши две функции, выразив x2 через x1 ---

# 1. Из x1^2/4 + x2^2 = 1  => x2 = sqrt(1 - x1^2/4)
# (Берем положительный корень, т.к. ищем решение в первой четверти)
def x2_from_ellipse(x1):
    # np.sqrt кидает ошибку на отрицательных числах, обрабатываем это
    # Вычисляем подкоренное выражение
    radicand = 1 - (x1**2) / 4.0
    # Возвращаем NaN (not a number) там, где корень не определен. Matplotlib не нарисует эти точки.
    return np.sqrt(np.where(radicand >= 0, radicand, np.nan))

# 2. Из 2*x2 - exp(x1) - x1 = 0  => x2 = (exp(x1) + x1) / 2
def x2_from_exp(x1):
    return (np.exp(x1) + x1) / 2.0

# --- Построение графика ---

# Создаем массив точек x1 от 0 до 2 (т.к. эллипс определен только на [-2, 2])
x1_vals = np.linspace(0, 2, 400)

# Вычисляем соответствующие значения x2 для каждой кривой
x2_vals_ellipse = x2_from_ellipse(x1_vals)
x2_vals_exp = x2_from_exp(x1_vals)

# Настройка и создание графика
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(10, 8))

# 1. Рисуем кривые
ax.plot(x1_vals, x2_vals_ellipse, label='$x_2 = \\sqrt{1 - x_1^2/4}$ (Эллипс)', color='blue', linewidth=2)
ax.plot(x1_vals, x2_vals_exp, label='$x_2 = (e^{x_1} + x_1)/2$', color='green', linewidth=2)

# 2. Отмечаем точку начального приближения
ax.plot(X1_INIT, X2_INIT, 'ro', markersize=9, label=f'Начальное приближение ({X1_INIT}, {X2_INIT})')
# Добавляем "крестик" в точку пересечения для наглядности
ax.axhline(X2_INIT, color='gray', linestyle=':', xmax=0.25)
ax.axvline(X1_INIT, color='gray', linestyle=':', ymax=0.92)

# 3. Оформление
ax.set_title('Графический поиск решения системы нелинейных уравнений', fontsize=16)
ax.set_xlabel('$x_1$', fontsize=14)
ax.set_ylabel('$x_2$', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True)
ax.axis('equal') # Делаем масштабы по осям одинаковыми, чтобы эллипс не искажался

# Устанавливаем пределы для лучшей видимости
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.2)

# Показать график
plt.show()