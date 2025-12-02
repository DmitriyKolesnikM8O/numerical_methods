import numpy as np
import matplotlib.pyplot as plt

# --- 1. ДАННЫЕ ИЗ ВЫВОДА ТВОЕЙ C++ ПРОГРАММЫ ---
# Я уже скопировал все необходимые значения из твоего вывода.

# Узлы сетки и точное решение
x_k = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
y_exact = np.array([
    1.0000000, 1.2825975, 1.5217091, 1.7068486, 1.8315322, 1.8937570,
    1.8960887, 1.8453513, 1.7519473, 1.6288697, 1.4904951
])

# Результаты y_k для каждого метода
y_results = {
    'euler_exp': [1.0, 1.3, 1.57, 1.7955819, 1.9649183, 2.0699109, 2.1070778, 2.0780866, 1.9898498, 1.8541346, 1.6866714],
    'euler_imp': [1.0, 1.2571265, 1.4620812, 1.6096196, 1.6986681, 1.7321007, 1.7162197, 1.6600078, 1.5742324, 1.4704890, 1.3602690],
    'euler_cauchy': [1.0, 1.285, 1.5258793, 1.7118925, 1.8364333, 1.8975336, 1.8979479, 1.8448161, 1.7489341, 1.6237021, 1.4838448],
    'euler_cauchy_iter': [1.0, 1.2815682, 1.5202008, 1.7055449, 1.8311563, 1.8949721, 1.8994063, 1.8510597, 1.7600711, 1.6391683, 1.5025002],
    'euler_imp_mid': [1.0, 1.285, 1.5257374, 1.7114171, 1.8354021, 1.8957159, 1.8951288, 1.8408192, 1.7436419, 1.6170709, 1.4759138],
    'rk3': [1.0, 1.2825330, 1.5215290, 1.7065234, 1.8310595, 1.8931628, 1.8954245, 1.8446865, 1.7513584, 1.6284281, 1.4902543],
    'rk4': [1.0, 1.2825928, 1.5217012, 1.7068393, 1.8315236, 1.8937507, 1.8960855, 1.8453511, 1.7519487, 1.6288700, 1.4904910],
    'adams_exp': [1.0, 1.2825928, 1.5217012, 1.7068393, 1.8312862, 1.8936834, 1.8963257, 1.8461092, 1.7534303, 1.6310913, 1.4934525],
    'adams_pc': [1.0, 1.2825928, 1.5217012, 1.7068393, 1.8315804, 1.8938441, 1.8961910, 1.8454399, 1.7519937, 1.6288500, 1.4903939]
}

# Имена методов для легенд
method_names = {
    'euler_exp': 'Явный Эйлер', 'euler_imp': 'Неявный Эйлер', 'euler_cauchy': 'Эйлер-Коши (явный)',
    'euler_cauchy_iter': 'Эйлер-Коши (неявный)', 'euler_imp_mid': 'Улучшенный Эйлер', 'rk3': 'Рунге-Кутта 3',
    'rk4': 'Рунге-Кутта 4', 'adams_exp': 'Адамс (явный)', 'adams_pc': 'Адамс (Предиктор-Корректор)'
}

# Данные для Рунге-Ромберга
final_errors = {
    'Явный Эйлер': {'rr': 0.08871508, 'true': 0.07161708},
    'Неявный Эйлер': {'rr': 0.08309249, 'true': 0.13022522},
    'Эйлер-Коши (явный)': {'rr': 0.02091312, 'true': 0.0066503},
    'Эйлер-Коши (неявный)': {'rr': 0.01173611, 'true': 0.01200508},
    'Улучшенный Эйлер': {'rr': 0.01323652, 'true': 0.01458132},
    'Рунге-Кутта 3': {'rr': 0.00038536, 'true': 0.0002408},
    'Рунге-Кутта 4': {'rr': 0.00000468, 'true': 0.0000041},
    'Адамс (явный)': {'rr': 0.00478980, 'true': 0.0029575},
    'Адамс (П-К)': {'rr': 0.00016711, 'true': 0.0001012}
}

# --- 2. НАСТРОЙКА СТИЛЯ ГРАФИКОВ ---
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 8) # Размер по умолчанию
plt.rcParams['font.size'] = 12           # Размер шрифта по умолчанию

# --- 3. ПОСТРОЕНИЕ ГРАФИКОВ ---

# == ГРАФИК 1: Погрешности методов низкого порядка ==
fig1, ax1 = plt.subplots()
ax1.axhline(0, color='black', linestyle='-', linewidth=1.5, label='Ноль (идеальное решение)') # Линия нуля для ориентира
low_order_methods = ['euler_exp', 'euler_imp', 'euler_cauchy', 'euler_cauchy_iter', 'euler_imp_mid']
for method in low_order_methods:
    # Вычисляем погрешность для каждого метода
    error = np.array(y_results[method]) - y_exact
    ax1.plot(x_k, error, 'o--', markersize=4, label=method_names[method])
ax1.set_title('Графики погрешностей: Методы низкого порядка (1 и 2)', fontsize=16)
ax1.set_xlabel('x', fontsize=12)
ax1.set_ylabel('Погрешность (y_k - y_ист)', fontsize=12)
ax1.legend()
ax1.grid(True)

# == ГРАФИК 2: Погрешности методов высокого порядка ==
fig2, ax2 = plt.subplots()
ax2.axhline(0, color='black', linestyle='-', linewidth=1.5, label='Ноль (идеальное решение)')
high_order_methods = ['rk3', 'rk4', 'adams_exp', 'adams_pc']
for method in high_order_methods:
    # Вычисляем погрешность для каждого метода
    error = np.array(y_results[method]) - y_exact
    ax2.plot(x_k, error, 'o--', markersize=5, label=method_names[method])
ax2.set_title('Графики погрешностей: Методы высокого порядка (3 и 4)', fontsize=16)
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel('Погрешность (y_k - y_ист)', fontsize=12)
ax2.legend()
ax2.grid(True)

# == ГРАФИК 3: Рост абсолютной ошибки (логарифмическая шкала) - ОСТАЕТСЯ БЕЗ ИЗМЕНЕНИЙ ==
# Этот график уже показывает то, что нужно, и делает это очень наглядно.
fig3, ax3 = plt.subplots()
for method_key, method_name in method_names.items():
    error = np.abs(np.array(y_results[method_key]) - y_exact)
    ax3.plot(x_k, error, 'o--', markersize=4, label=method_name)
ax3.set_yscale('log')
ax3.set_title('Рост абсолютной погрешности |y_ист - y_k|', fontsize=16)
ax3.set_xlabel('x', fontsize=12)
ax3.set_ylabel('Погрешность (log шкала)', fontsize=12)
ax3.legend(loc='upper left', fontsize=10)
ax3.grid(True, which='both')


# == ГРАФИК 4: Сравнение оценки по Рунге-Ромбергу с истинной ошибкой ==
labels = list(final_errors.keys())
rr_errors = [d['rr'] for d in final_errors.values()]
true_errors = [d['true'] for d in final_errors.values()]

x_pos = np.arange(len(labels))
width = 0.35

fig4, ax4 = plt.subplots(figsize=(14, 9))
rects1 = ax4.bar(x_pos - width/2, rr_errors, width, label='Оценка по Рунге-Ромбергу')
rects2 = ax4.bar(x_pos + width/2, true_errors, width, label='Истинная ошибка')

ax4.set_yscale('log')
ax4.set_ylabel('Погрешность в конечной точке (log шкала)', fontsize=12)
ax4.set_title('Сравнение оценки по Рунге-Ромбергу с истинной ошибкой', fontsize=16)
ax4.set_xticks(x_pos)
ax4.set_xticklabels(labels, rotation=45, ha="right")
ax4.legend()
ax4.grid(True, which='both', axis='y')


# == НОВЫЙ ГРАФИК 5: "Под микроскопом" - сравнение точных методов в конце отрезка ==
fig5, ax5 = plt.subplots()

# Строим точное решение, но делаем его более гладким (больше точек)
x_smooth = np.linspace(x_k[-3], x_k[-1], 300) # Берем последние 3 точки, т.е. [0.8, 1.0]
y_smooth_exact = np.cos(x_smooth)**3 + np.sin(x_smooth)*(1 + 2*np.cos(x_smooth)**2)
ax5.plot(x_smooth, y_smooth_exact, 'k-', linewidth=4, alpha=0.8, label='Точное решение')

# Выбираем только самые точные методы для сравнения
most_accurate_methods = ['rk4', 'adams_pc']
for method in most_accurate_methods:
    # Берем только последние 3 точки из результатов
    ax5.plot(x_k[-3:], y_results[method][-3:], 'o--', markersize=6, label=method_names[method])

ax5.set_title('Сравнение точных методов на отрезке [0.8, 1.0]', fontsize=16)
ax5.set_xlabel('x', fontsize=12)
ax5.set_ylabel('y(x)', fontsize=12)
ax5.legend()
ax5.grid(True)


# Показываем все созданные графики
plt.tight_layout()
plt.show()