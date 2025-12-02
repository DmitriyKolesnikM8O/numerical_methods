import numpy as np
import matplotlib.pyplot as plt

# --- 1. ДАННЫЕ ИЗ ВЫВОДА ТВОЕЙ C++ ПРОГРАММЫ (перепроверенные) ---

# Узлы сетки
x_k = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

# Точное решение
y_exact = np.array([
    0.50000000, 0.58019173, 0.67156876, 0.77541629, 0.89313704,
    1.02626194, 1.17646249, 1.34556494, 1.53556645, 1.74865361, 1.98722325
])

# Численное решение: Метод стрельбы
y_shooting = np.array([
    0.50000647, 0.58019804, 0.67157490, 0.77542225, 0.89314281,
    1.02626752, 1.17646787, 1.34557009, 1.53557137, 1.74865828, 1.98722764
])

# Численное решение: Конечно-разностный метод
y_finite_diff = np.array([
    0.49905088, 0.57904850, 0.67021675, 0.77383888, 0.89131528,
    1.02417413, 1.17408372, 1.34286651, 1.53251537, 1.74521193, 1.98334739
])

# Данные для Рунге-Ромберга и истинной ошибки
# ИСПРАВЛЕНЫ ЗНАЧЕНИЯ
final_errors = {
    'Метод стрельбы': {'rr': 0.00000580, 'true_max': 0.00000647},
    'Конечно-разностный метод': {'rr': 0.00384302, 'true_max': 0.00387586}
}

# --- 2. НАСТРОЙКА СТИЛЯ ГРАФИКОВ ---
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 12

# --- 3. ПОСТРОЕНИЕ ГРАФИКОВ ---

# == ГРАФИК 1: Сравнение численных решений с точным ==
fig1, ax1 = plt.subplots(figsize=(12, 8))
ax1.plot(x_k, y_exact, 'k-', linewidth=4, alpha=0.8, label='Точное решение')
ax1.plot(x_k, y_shooting, 'o--', color='blue', markersize=6, label='Метод стрельбы')
ax1.plot(x_k, y_finite_diff, 's:', color='red', markersize=5, label='Конечно-разностный метод')
ax1.set_title('Сравнение численных решений краевой задачи', fontsize=16)
ax1.set_xlabel('x', fontsize=12)
ax1.set_ylabel('y(x)', fontsize=12)
ax1.legend()
ax1.grid(True, which='both')

# == ГРАФИК 2: Графики абсолютных погрешностей ==
error_shooting = np.abs(y_shooting - y_exact)
error_finite_diff = np.abs(y_finite_diff - y_exact)

fig2, ax2 = plt.subplots(figsize=(12, 8))
ax2.plot(x_k, error_shooting, 'o-', color='blue', label=f'Метод стрельбы (max_err={np.max(error_shooting):.2e})')
ax2.plot(x_k, error_finite_diff, 's-', color='red', label=f'Конечно-разностный метод (max_err={np.max(error_finite_diff):.2e})')
ax2.set_title('Абсолютная погрешность |y_ист - y_числ|', fontsize=16)
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel('Погрешность', fontsize=12)
# Используем логарифмическую шкалу, чтобы увидеть обе линии
ax2.set_yscale('log')
ax2.legend()
ax2.grid(True, which='both')

# == ГРАФИК 3: Сравнение оценки по Рунге-Ромбергу с истинной ошибкой ==
labels = list(final_errors.keys())
rr_errors = [d['rr'] for d in final_errors.values()]
true_errors = [d['true_max'] for d in final_errors.values()]

x_pos = np.arange(len(labels))
width = 0.35

fig3, ax3 = plt.subplots(figsize=(10, 7))
rects1 = ax3.bar(x_pos - width/2, rr_errors, width, label='Оценка по Рунге-Ромбергу')
rects2 = ax3.bar(x_pos + width/2, true_errors, width, label='Максимальная истинная ошибка')

ax3.set_ylabel('Погрешность', fontsize=12)
ax3.set_title('Сравнение оценки по Рунге-Ромбергу с истинной ошибкой', fontsize=16)
ax3.set_xticks(x_pos)
ax3.set_xticklabels(labels)
ax3.legend()
# Логарифмическая шкала здесь тоже поможет увидеть разницу в масштабах
ax3.set_yscale('log')
ax3.grid(True, which='both', axis='y')

# Показываем все созданные графики
plt.tight_layout()
plt.show()