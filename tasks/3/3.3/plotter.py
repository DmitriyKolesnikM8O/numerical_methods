import numpy as np
import matplotlib.pyplot as plt



x_data = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
y_data = [-0.4597, 1.0, 1.5403, 1.5839, 2.010, 3.3464]
a0_1 = 0.55616 
a1_1 = 0.63155
a0_2 = 0.56894
a1_2 = 0.68904
a2_2 = -0.01917

# методичка
# x_data = [0.0, 1.7, 3.4, 5.1, 6.8, 8.5]
# y_data = [0.0, 1.30380, 1.84390, 2.25830, 2.60770, 2.91550]
# a0_1 = 0.47128 
# a1_1 = 0.31771
# a0_2 = 0.12944
# a1_2 = 0.61933
# a2_2 = -0.03548


# другой вариант
# x_data = [-0.9, 0.0, 0.9, 1.8, 2.7, 3.6]
# y_data = [-0.36892, 0.0, 0.36892, 0.85408, 1.78560, 6.31380]
# a0_1 = -0.19013
# a1_1 = 1.24621
# a0_2 = -0.46450
# a1_2 = -0.12563
# a2_2 = 0.50809

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