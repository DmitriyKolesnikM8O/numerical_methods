import pandas as pd
import matplotlib.pyplot as plt

# Читаем данные из CSV
data = pd.read_csv('./tasks/3/3.1/plot_data.csv')

# Строим графики
plt.figure(figsize=(10, 6))
plt.plot(data['x'], data['actual_f(x)'], label='Исходная функция f(x) = cos(x) + x', linewidth=2)
plt.plot(data['x'], data['polynomial_p(x)'], label='Интерполяционный многочлен P(x)', linestyle='--')

# Отмечаем узлы интерполяции на графике
nodes_x = [0, 3.14159/6, 3.14159/4, 3.14159/2]
nodes_y = [1.0, 1.389711, 1.492474, 1.570796]
plt.plot(nodes_x, nodes_y, 'ro', label='Узлы интерполяции')

plt.title('Сравнение функции и интерполяционного многочлена')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show()
