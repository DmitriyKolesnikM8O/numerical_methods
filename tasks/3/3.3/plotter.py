import numpy as np
import matplotlib.pyplot as plt


# 1. Исходные табличные данные (Вариант 13)
x_data = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
y_data = [-0.4597, 1.0, 1.5403, 1.5839, 2.010, 3.3464]

# 2. Коэффициенты для многочлена 1-й степени (F1(x) = a0 + a1*x)
a0_1 = 0.55616 
a1_1 = 0.63155

# 3. Коэффициенты для многочлена 2-й степени (F2(x) = a0 + a1*x + a2*x^2)
a0_2 = 0.56894
a1_2 = 0.68904
a2_2 = -0.01917

def p1(x):
    return a0_1 + a1_1 * x

def p2(x):
    return a0_2 + a1_2 * x + a2_2 * x**2

x_range = np.linspace(min(x_data), max(x_data), 200)
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(x_data, y_data, 'o', color='black', label='табличные данные')
ax.plot(x_range, p1(x_range), color='black', linestyle='-', label='приближающий многочлен первой степени')
ax.plot(x_range, p2(x_range), color='black', linestyle='--', label='приближающий многочлен второй степени')
ax.set_title('Графики приближающих многочленов (Вариант 13)', fontsize=14)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.legend()
ax.grid(True)


plt.show()