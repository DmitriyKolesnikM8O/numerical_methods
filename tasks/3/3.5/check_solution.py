import numpy as np
from scipy.integrate import quad

def f(x):
    return (x**2) / (x**3 - 27)

a = -2.0
b = 2.0


result, error_estimate = quad(f, a, b)

print("--- Вычисление точного значения интеграла ---")
print(f"Интеграл от f(x) = x^2 / (x^3 - 27) на отрезке [{a}, {b}]")
print("-" * 50)
print(f"Точное значение: {result:.7f}")
print(f"Оценка погрешности вычисления (SciPy): {error_estimate}")