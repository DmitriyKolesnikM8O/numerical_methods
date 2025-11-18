import numpy as np
from scipy.interpolate import CubicSpline

x_nodes = [0.0, 1.0, 2.0, 3.0, 4.0]
y_nodes = [1.0, 1.5403, 1.5839, 2.010, 3.3464]

spline = CubicSpline(x_nodes, y_nodes, bc_type='natural')

x_star = 1.5

y_star = spline(x_star)

print(f"Точки X: {x_nodes}")
print(f"Точки Y: {y_nodes}")
print(f"Граничные условия: 'natural' (вторая производная на концах равна 0)")
print("-" * 40)
print(f"Значение сплайна в точке x = {x_star}:")
print(f"S({x_star}) = {y_star:.8f}")