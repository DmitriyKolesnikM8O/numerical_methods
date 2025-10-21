import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.log(x + 1) - 2 * x + 0.5

x_vals = np.linspace(-0.5, 2, 400) # Широкий диапазон
y_vals = f(x_vals)

plt.plot(x_vals, y_vals)
plt.grid(True)
plt.axhline(0, color='black')
plt.show()