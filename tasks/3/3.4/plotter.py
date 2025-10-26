import numpy as np
import matplotlib.pyplot as plt

x_nodes = [0.2, 0.5, 0.8, 1.1, 1.4]
y_nodes = [12.906, 5.5273, 3.8777, 3.2692, 3.0319]

x_star = 0.8

y_prime_star = -3.76350  
y_double_prime_star = 11.56778 

idx_star = x_nodes.index(x_star)
x_central = x_nodes[idx_star-1 : idx_star+2]
y_central = y_nodes[idx_star-1 : idx_star+2]

coeffs = np.polyfit(x_central, y_central, 2)
a, b, c = coeffs

def parabola(x):
    return a * x**2 + b * x + c

def first_derivative(x):
    return 2 * a * x + b

def second_derivative(x):

    return 2 * a

x_range = np.linspace(x_central[0], x_central[2], 200)

tangent_line = y_prime_star * (x_range - x_star) + parabola(x_star)

plt.style.use('seaborn-v0_8-whitegrid')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
fig.suptitle('Визуализация численного дифференцирования в точке x* = {}'.format(x_star), fontsize=16)
ax1.set_title('Аппроксимация функции параболой', fontsize=14)
ax1.plot(x_nodes, y_nodes, 'o', color='gray', markersize=8, label='Исходные точки')
ax1.plot(x_central, y_central, 'o', color='red', markersize=10, label='Узлы для аппроксимации')
ax1.plot(x_range, parabola(x_range), color='blue', linewidth=2, label='Аппроксимирующая парабола')
ax1.plot(x_range, tangent_line, '--', color='green', linewidth=2, label=f'Касательная (y\' = {y_prime_star:.4f})')
ax1.set_ylabel('y(x)', fontsize=12)
ax1.legend()
ax1.grid(True)

ax2.set_title('Графики производных на отрезке', fontsize=14)
ax2.plot(x_range, first_derivative(x_range), color='green', linewidth=2, label=f'Первая производная y\'(x)')
ax2.plot(x_range, [second_derivative(x) for x in x_range], '--', color='purple', linewidth=2, label=f'Вторая производная y\'\'(x) = {y_double_prime_star:.4f}')
ax2.plot(x_star, y_prime_star, 'go')
ax2.plot(x_star, y_double_prime_star, 'o', color='purple')
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel("y'(x), y''(x)", fontsize=12)
ax2.legend()
ax2.grid(True)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()