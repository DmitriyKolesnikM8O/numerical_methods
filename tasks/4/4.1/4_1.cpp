#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <functional>
#include <string>
#include "Vector.hpp"

// Определяем типы для наших функций
using ODEFunc = std::function<double(double, double, double)>;
using ExactSolutionFunc = std::function<double(double)>;

class Task_4_1 {
public:
    Task_4_1(ODEFunc f_y, ODEFunc f_z, ExactSolutionFunc y_exact_func,
               double x_start, double y_start, double z_start,
               double x_end, double step)
        : fy(f_y), fz(f_z), y_exact(y_exact_func),
          x0(x_start), y0(y_start), z0(z_start),
          xk(x_end), h(step) {}

    void Do() {
        std::cout << std::fixed << std::setprecision(7);
        PrintInitialData();
        
        // Методы, которые выводят пошаговые таблицы, остаются без изменений
        SolveEulerExplicit();
        SolveEulerImplicit();
        SolveEulerCauchy();
        SolveEulerCauchyIterative();
        SolveEulerImprovedMidpoint();
        SolveRK3();
        SolveRK4();
        SolveAdamsExplicit();
        SolveAdamsPredictorCorrector();

        // А эта функция теперь будет работать для всех методов
        EvaluateRungeRombergForAll();
    }

private:
    ODEFunc fy, fz;
    ExactSolutionFunc y_exact;
    double x0, y0, z0, xk, h;

    // =========================================================================
    // ОСНОВНЫЕ МЕТОДЫ РЕШЕНИЯ (с выводом таблиц) - БЕЗ ИЗМЕНЕНИЙ
    // =========================================================================
    
    // 1. Явный метод Эйлера (с печатью)
    void SolveEulerExplicit() {
        PrintTableHeader("1. Явный метод Эйлера");
        double x = x0, y = y0, z = z0; int n = static_cast<int>(round((xk-x0)/h));
        for (int k = 0; k <= n; ++k) {
            if (k == n) { PrintTableRow(k,x,y,z); break; }

            // 1. Вычисляем приращения dy_k и dz_k.
            //    Это применение общей формулы Эйлера (4.2): y_{k+1} = y_k + h*f(x_k, y_k)
            //    Мы применяем ее для каждой из функций нашей системы: y' = z и z' = f(x,y,z).

            // Приращение для y: dy_k = h * y'_k = h * z_k
            double dy = h * fy(x,y,z); // fy(x,y,z) возвращает z

            // Приращение для z: dz_k = h * z'_k = h * f(x_k, y_k, z_k)
            double dz = h * fz(x,y,z); // fz(...) возвращает 2*tg(x)*z + 3*y
            PrintTableRow(k, x, y, z, dy, dz);

            // 3. Вычисляем значения на следующем шаге k+1.
            //    y_{k+1} = y_k + dy_k
            //    z_{k+1} = z_k + dz_k
            //    x_{k+1} = x_k + h
            y += dy; z += dz; x += h;
        }
    }
    // 2. Неявный метод Эйлера (с печатью)
    void SolveEulerImplicit() {
        PrintTableHeader("2. Неявный метод Эйлера (решен итерациями)");
        double x = x0, y = y0, z = z0; int n = static_cast<int>(round((xk-x0)/h));
        for (int k = 0; k <= n; ++k) {
            if (k == n) { PrintTableRow(k,x,y,z); break; }
            
            double next_x = x + h;
            // Прогноз
            // ЦЕЛЬ: Решить систему неявных уравнений Эйлера (формула 4.3):
            //   y_{k+1} = y_k + h * f_y(x_{k+1}, y_{k+1}, z_{k+1})
            //   z_{k+1} = z_k + h * f_z(x_{k+1}, y_{k+1}, z_{k+1})
            // где f_y = z, а f_z = 2*tg(x)*z + 3y.
            //
            // МЕТОД РЕШЕНИЯ: Метод простой итерации (как описано на стр. 4 под формулой 4.6).

            // 1. "Прогноз" (начальное приближение для итераций).
            //    Чтобы начать итерационный процесс, нам нужно первое предположение
            //    о том, чему могут быть равны y_{k+1} и z_{k+1}.
            //    Самый простой способ - сделать шаг явным методом Эйлера.
            double y_next = y + h*fy(x,y,z);
            double z_next = z + h*fz(x,y,z);
            // Итерации для поиска y_{k+1} и z_{k+1}

            // 2. Итерационное уточнение.
            //    Мы несколько раз (здесь - 10) подставляем текущее приближение
            //    в правую часть неявного уравнения, чтобы получить новое, более точное.
            //    Это реализация метода Эйлера-Коши с итерационной обработкой,
            //    примененная для решения более простого неявного уравнения Эйлера.
            for(int i=0;i<10;++i){

                // ВАЖНО: для решения системы уравнений f_y(..., z_{k+1}) и f_z(..., y_{k+1}, ...)
                // итерации лучше делать поочередно (метод Якоби или Зейделя).
                // Здесь используется вариант, похожий на метод Зейделя, где новое значение
                // одной переменной сразу используется для вычисления другой.
                // Чтобы избежать расходимости, для f_z используем y_next от предыдущей итерации.
                double y_iter_old = y_next;
                // y_new = y_k + h * z_old
                y_next = y + h * fy(next_x, y_next, z_next);

                // z_new = z_k + h * f_z(x_{k+1}, y_old_from_this_step, z_old)
                z_next = z + h * fz(next_x, y_iter_old, z_next);
            }

            // ВЫЧИСЛЯЕМ ПРИРАЩЕНИЯ ПОСТФАКТУМ
            // 3. Вычисляем итоговые приращения.
            //    dy_k - это фактическое изменение 'y' за весь шаг.
            double dy = y_next - y;
            double dz = z_next - z;
            PrintTableRow(k, x, y, z, dy, dz); // Печатаем с приращениями
            
            // 5. Обновляем переменные до состояния k+1.
            y = y_next; z = z_next; x = next_x;
        }
    }
    // 3. Метод Эйлера-Коши (с печатью)
    void SolveEulerCauchy() {
        PrintTableHeader("3. Метод Эйлера-Коши (явный Предиктор-Корректор)");
        double x = x0, y = y0, z = z0; int n = static_cast<int>(round((xk-x0)/h));
        for (int k = 0; k <= n; ++k) {
            if (k == n) { PrintTableRow(k,x,y,z); break; }
            double next_x = x + h;

            // 1. ЭТАП ПРЕДИКТОРА ("прогноз").
            //    Сначала делаем грубый шаг обычным методом Эйлера, чтобы получить
            //    предварительное значение в конце интервала (y_pred, z_pred).
            //    Это соответствует `ỹ_{k+1}` из формулы (4.4) в методичке.
            double y_pred = y + h * fy(x,y,z), z_pred = z + h * fz(x,y,z);

            // 2. ЭТАП КОРРЕКТОРА ("уточнение").
            //    Вычисляем финальные приращения dy и dz, используя ПОЛУСУММУ (среднее арифметическое)
            //    наклонов в НАЧАЛЕ интервала (f(x_k, y_k, z_k)) и в ПРЕДСКАЗАННОМ КОНЦЕ
            //    интервала (f(x_{k+1}, y_pred, z_pred)).
            //    Это в точности реализует правую часть формулы (4.4).
            
            // Приращение для y: dy_k = (h/2) * [y'_k + y'_{k+1, pred}] = (h/2) * [z_k + z_pred]
            // f_y(x,y,z)=z, f_y(next_x,y_pred,z_pred)=z_pred
            double dy = (h/2.0) * (fy(x,y,z) + fy(next_x, y_pred, z_pred));

            // Приращение для z: dz_k = (h/2) * [z'_k + z'_{k+1, pred}]
            double dz = (h/2.0) * (fz(x,y,z) + fz(next_x, y_pred, z_pred));
            PrintTableRow(k, x, y, z, dy, dz);
            y += dy; z += dz; x = next_x;
        }
    }
    // 4. Неявный Эйлер-Коши (с печатью)
    void SolveEulerCauchyIterative() {
        PrintTableHeader("4, 5. Метод Эйлера-Коши с итерационной обработкой");
        double x = x0, y = y0, z = z0; int n = static_cast<int>(round((xk-x0)/h));
        for (int k = 0; k <= n; ++k) {
            if (k == n) { PrintTableRow(k,x,y,z); break; }

            double next_x = x + h;
            // Прогноз

            // ---= ШАГ ВЫЧИСЛЕНИЙ СОГЛАСНО МЕТОДИЧКЕ (формула 4.6) =---

            // ЦЕЛЬ: Решить систему неявных уравнений Эйлера-Коши (метод трапеций, формула 4.5):
            //   y_{k+1} = y_k + (h/2) * [ f_y(x_k, y_k, z_k) + f_y(x_{k+1}, y_{k+1}, z_{k+1}) ]
            //   z_{k+1} = z_k + (h/2) * [ f_z(x_k, y_k, z_k) + f_z(x_{k+1}, y_{k+1}, z_{k+1}) ]
            //
            // МЕТОД РЕШЕНИЯ: Метод простой итерации, как в формуле (4.6).

            // 1. "Прогноз" (начальное приближение y_{k+1}^{(0)} для итераций).
            //    Как и в неявном методе Эйлера, нам нужно первое предположение.
            //    Используем явный метод Эйлера.
            double y_next = y + h*fy(x,y,z);
            double z_next = z + h*fz(x,y,z);

            // 2. Итерационное уточнение (реализация формулы 4.6).
            //    Заранее вычисляем наклоны в НАЧАЛЬНОЙ точке, т.к. они не меняются в цикле.
            double fy0 = fy(x,y,z), fz0 = fz(x,y,z);

            // Итерации для поиска y_{k+1} и z_{k+1}
            // Запускаем цикл для "подгонки" y_{k+1} и z_{k+1}.
            for(int i=0; i<10; ++i){
                // Сохраняем y_next от предыдущей итерации для стабильности решения системы
                double y_iter_old = y_next;

                // y_{k+1}^{(i)} = y_k + (h/2) * [ f_y(x_k,...) + f_y(x_{k+1}, y_{k+1}^{(i-1)}, ...) ]
                y_next = y + (h/2.0) * (fy0 + fy(next_x, y_next, z_next));

                // z_{k+1}^{(i)} = z_k + (h/2) * [ f_z(x_k,...) + f_z(x_{k+1}, y_{k+1}^{(i-1)}, ...) ]
                z_next = z + (h/2.0) * (fz0 + fz(next_x, y_iter_old, z_next));
            }

            // ВЫЧИСЛЯЕМ ПРИРАЩЕНИЯ ПОСТФАКТУМ
            double dy = y_next - y;
            double dz = z_next - z;
            PrintTableRow(k, x, y, z, dy, dz); // Печатаем с приращениями
            
            y = y_next; z = z_next; x = next_x;
        }
    }
    // 6. Улучшенный метод Эйлера (с печатью)
   void SolveEulerImprovedMidpoint() {
        PrintTableHeader("6. Первый улучшенный метод Эйлера");
        double x = x0, y = y0, z = z0; int n = static_cast<int>(round((xk-x0)/h));
        for (int k = 0; k <= n; ++k) {
            if (k == n) { PrintTableRow(k,x,y,z); break; }

            // ---= ШАГ ВЫЧИСЛЕНИЙ СОГЛАСНО МЕТОДИЧКЕ (формула 4.7) =---
            // Идея: вместо усреднения наклонов по краям, как в методе Эйлера-Коши,
            // мы делаем шаг, используя более точный наклон, вычисленный в СЕРЕДИНЕ интервала.

            // 1. Находим координаты "средней точки" (x_{k+1/2}, y_{k+1/2}, z_{k+1/2}).
            //    Это делается с помощью "полушага" обычным методом Эйлера.
            double mid_x = x + h/2.0;

            // y_{k+1/2} = y_k + (h/2) * y'_k
            double y_mid = y + (h/2.0) * fy(x,y,z);

            // z_{k+1/2} = z_k + (h/2) * z'_k
            double z_mid = z + (h/2.0) * fz(x,y,z);

            // 2. Вычисляем финальные приращения dy и dz для ПОЛНОГО шага h.
            //    В качестве наклона используется производная, вычисленная в найденной
            //    "средней точке". Это вторая строка из формулы (4.7).
            
            // Приращение для y: dy_k = h * y'_{k+1/2} = h * z_{k+1/2}
            double dy = h * fy(mid_x, y_mid, z_mid);

            // Приращение для z: dz_k = h * z'_{k+1/2} = h * f(x_{k+1/2}, y_{k+1/2}, z_{k+1/2})
            double dz = h * fz(mid_x, y_mid, z_mid);
            PrintTableRow(k, x, y, z, dy, dz);
            y += dy; z += dz; x += h;
        }
    }
    // 7. Метод Рунге-Кутты 3-го порядка (с печатью)
    void SolveRK3() {
        PrintTableHeader("7. Метод Рунге-Кутты 3-го порядка");
        double x = x0, y = y0, z = z0; int n = static_cast<int>(round((xk-x0)/h));
        for (int k=0; k<=n; ++k) {
            if (k==n) { PrintTableRow(k,x,y,z); break; }

            // ---= ШАГ ВЫЧИСЛЕНИЙ СОГЛАСНО МЕТОДИЧКЕ (формула 4.9) =---
            // Идея: Вместо двух "пробных" наклонов, как в методах 2-го порядка,
            // здесь вычисляются три наклона (K1, K2, K3), чтобы получить
            // более точную оценку среднего наклона на интервале.
            // Для системы уравнений мы параллельно считаем K_i для 'y' и L_i для 'z'.

            // 1. Вычисляем наклоны в начальной точке (k_1).
            //    K_1 = h * y'(x_k, y_k, z_k)
            //    L_1 = h * z'(x_k, y_k, z_k)
            double K1=h*fy(x,y,z), L1=h*fz(x,y,z);

            // 2. Вычисляем наклоны в точке (x_k + h/3), используя K_1, L_1 (k_2).
            //    K_2 = h * y'(x_k + h/3, y_k + K_1/3, z_k + L_1/3)
            //    L_2 = h * z'(x_k + h/3, y_k + K_1/3, z_k + L_1/3)
            double K2=h*fy(x+h/3.,y+K1/3.,z+L1/3.), L2=h*fz(x+h/3.,y+K1/3.,z+L1/3.);

            // 3. Вычисляем наклоны в точке (x_k + 2h/3), используя K_2, L_2 (k_3).
            //    K_3 = h * y'(x_k + 2h/3, y_k + 2*K_2/3, z_k + 2*L_2/3)
            //    L_3 = h * z'(x_k + 2h/3, y_k + 2*K_2/3, z_k + 2*L_2/3)
            double K3=h*fy(x+2.*h/3.,y+2.*K2/3.,z+2.*L2/3.), L3=h*fz(x+2.*h/3.,y+2.*K2/3.,z+2.*L2/3.);

            // 4. Вычисляем итоговые приращения dy и dz как взвешенную сумму наклонов.
            //    Это в точности соответствует формуле (4.9): Δy_k = (1/4)*(K_1 + 3*K_3)
            double dy = (K1 + 3.*K3)/4.;
            double dz = (L1 + 3.*L3)/4.;
            PrintTableRow(k, x, y, z, dy, dz);
            y += dy; z += dz; x += h;
        }
    }
    // 8. Метод Рунге-Кутты 4-го порядка (с печатью)
    void SolveRK4() {
        PrintTableHeader("8. Метод Рунге-Кутты 4-го порядка");
        double x=x0, y=y0, z=z0; int n=static_cast<int>(round((xk-x0)/h));
        for (int k=0; k<=n; ++k) {
            if (k == n) { // Для последней строки просто печатаем результат
                PrintTableRow(k, x, y, z);
                break;
            }
            
            // Вычисляем коэффициенты
            // ---= ШАГ ВЫЧИСЛЕНИЙ СОГЛАСНО МЕТОДИЧКЕ (формула 4.10 и 4.16) =---
            // Идея: вычисляются четыре "пробных" наклона на интервале, которые
            // затем усредняются с весами для получения очень точного приращения.
            // Для системы уравнений мы параллельно считаем K_i для 'y' и L_i для 'z'.

            // 1. Вычисляем наклоны в НАЧАЛЬНОЙ точке (k_1).
            //    K_1 = h * y'(x_k, y_k, z_k)
            //    L_1 = h * z'(x_k, y_k, z_k)
            double K1=h*fy(x,y,z), L1=h*fz(x,y,z);

            // 2. Вычисляем наклоны в СЕРЕДИНЕ интервала, используя K_1, L_1 (k_2).
            //    K_2 = h * y'(x_k + h/2, y_k + K_1/2, z_k + L_1/2)
            //    L_2 = h * z'(x_k + h/2, y_k + K_1/2, z_k + L_1/2)
            double K2=h*fy(x+h/2,y+K1/2,z+L1/2), L2=h*fz(x+h/2,y+K1/2,z+L1/2);

            // 3. Вычисляем УТОЧНЕННЫЕ наклоны в СЕРЕДИНЕ интервала, используя K_2, L_2 (k_3).
            //    K_3 = h * y'(x_k + h/2, y_k + K_2/2, z_k + L_2/2)
            //    L_3 = h * z'(x_k + h/2, y_k + K_2/2, z_k + L_2/2)
            double K3=h*fy(x+h/2,y+K2/2,z+L2/2), L3=h*fz(x+h/2,y+K2/2,z+L2/2);

            // 4. Вычисляем наклоны в КОНЦЕ интервала, используя K_3, L_3 (k_4).
            //    K_4 = h * y'(x_k + h, y_k + K_3, z_k + L_3)
            //    L_4 = h * z'(x_k + h, y_k + K_3, z_k + L_3)
            double K4=h*fy(x+h,y+K3,z+L3), L4=h*fz(x+h,y+K3,z+L3);

            // --- НОВЫЙ БЛОК: Вычисление theta по формуле (4.17) ---
            double theta1 = 0.0, theta2 = 0.0;
            // Проверка на деление на ноль
            if (std::abs(K1 - K2) > 1e-12) {
                theta1 = std::abs((K2 - K3) / (K1 - K2));
            }
            if (std::abs(L1 - L2) > 1e-12) {
                theta2 = std::abs((L2 - L3) / (L1 - L2));
            }
            // Здесь можно было бы добавить вывод theta в таблицу
            // и логику изменения шага h
            std::cout << "Theta1 = " << theta1 << ", Theta2 = " << theta2 << "\n";
            
            // Вычисляем ПРИРАЩЕНИЯ dy и dz
            // 5. Вычисляем итоговые приращения dy и dz как взвешенную сумму наклонов.
            //    Это в точности соответствует формуле (4.10) (и 4.16 для системы):
            //    Δy_k = (1/6) * (K_1 + 2*K_2 + 2*K_3 + K_4)
            double dy = (K1 + 2*K2 + 2*K3 + K4) / 6.0;
            double dz = (L1 + 2*L2 + 2*L3 + L4) / 6.0;
            
            // Печатаем строку, включая dy и dz
            PrintTableRow(k, x, y, z, dy, dz);

            // Обновляем значения для следующего шага
            y += dy;
            z += dz;
            x += h;
        }
    }
    
    // 9. Метод Адамса (явный) (с печатью)
    void SolveAdamsExplicit() {
        // Заголовок для Адамса не включает dy/dz, т.к. приращение явно не вычисляется.
        // Вместо этого используется PrintTableRow без этих аргументов.
        PrintTableHeader("9. Метод Адамса (явный, 4-го порядка)"); // Используем PrintTableRow без dy, dz

        // Метод требует 4 начальные точки, т.е. минимум 4 интервала.
        int n=static_cast<int>(round((xk-x0)/h)); if(n<4){std::cout<<"Мало шагов\n"; return;}

        // ---= ЭТАП 1: "РАЗГОН" (Подготовка истории) =---
        // Метод Адамса - многошаговый, ему нужна "история" предыдущих точек.
        // Как указано в методичке (стр. 19), первые точки (y1, y2, y3)
        // необходимо получить с помощью одношагового метода, например, РК4.

        // Создаем векторы для хранения истории: x_k, y_k, z_k и их производных f_y, f_z
        Vector xh(n+1), yh(n+1), zh(n+1), fyh(n+1), fzh(n+1);

        // Вызываем вспомогательную функцию, которая заполняет первые 4 элемента
        // этих векторов (индексы 0, 1, 2, 3), используя РК4.
        RungeKuttaStartup(xh, yh, zh, fyh, fzh, h);

        for(int k=0; k<4; ++k) PrintTableRow(k, xh.Get(k), yh.Get(k), zh.Get(k));

        // ---= ЭТАП 2: ОСНОВНОЙ ЦИКЛ МЕТОДА АДАМСА =---
        // Начинаем с k=3, чтобы вычислить y_4, используя историю из точек 0, 1, 2, 3.
        for(int k=3; k<n; ++k) {
            // Получаем значения y_k и z_k из истории
            double current_y = yh.Get(k);
            double current_z = zh.Get(k);

            // Считаем y_{k+1} и z_{k+1}
            // Применяем явную формулу Адамса (формула 4.25 из методички)
            // для вычисления y_{k+1} и z_{k+1}.
            // y_{k+1} = y_k + (h/24) * (55*f_y,k - 59*f_y,k-1 + 37*f_y,k-2 - 9*f_y,k-3)
            double next_y = current_y + (h/24.0) * (55*fyh.Get(k) - 59*fyh.Get(k-1) + 37*fyh.Get(k-2) - 9*fyh.Get(k-3));

            // z_{k+1} = z_k + (h/24) * (55*f_z,k - 59*f_z,k-1 + 37*f_z,k-2 - 9*f_z,k-3)
            double next_z = current_z + (h/24.0) * (55*fzh.Get(k) - 59*fzh.Get(k-1) + 37*fzh.Get(k-2) - 9*fzh.Get(k-3));
            double next_x = xh.Get(k) + h;

            // Сохраняем в историю
            yh.Set(k+1, next_y);
            zh.Set(k+1, next_z);
            xh.Set(k+1, next_x);
            
            // ВАЖНО: Считаем ОБЕ производные от НОВЫХ, согласованных (y,z)
            fyh.Set(k+1, fy(next_x, next_y, next_z));
            fzh.Set(k+1, fz(next_x, next_y, next_z));
            
            PrintTableRow(k+1, next_x, next_y, next_z);
        }
    }

    // 10. Метод Адамса (предиктор-корректор)
    void SolveAdamsPredictorCorrector() {
        // Заголовок для Адамса не включает dy/dz. Вместо этого в методичке
        // (Таблица 4.14) показаны предсказанные и скорректированные значения.
        // Мы выводим только финальные (скорректированные).
        PrintTableHeader("10. Метод Адамса-Бэшфортса-Моултона");
        int n=static_cast<int>(round((xk-x0)/h)); if(n<4){std::cout<<"Мало шагов\n"; return;}

        // ---= ЭТАП 1: "РАЗГОН" (Подготовка истории) =---
        // Аналогично предыдущему методу, для старта нам нужны первые 4 точки (0,1,2,3),
        // вычисленные более точным одношаговым методом (РК4).
        Vector xh(n+1), yh(n+1), zh(n+1), fyh(n+1), fzh(n+1);
        RungeKuttaStartup(xh, yh, zh, fyh, fzh, h);

        // Печатаем "разогнанные" точки.
        for(int k=0; k<4; ++k) PrintTableRow(k, xh.Get(k), yh.Get(k), zh.Get(k));
        
        // ---= ЭТАП 2: ОСНОВНОЙ ЦИКЛ ПРЕДИКТОР-КОРРЕКТОР =---
        // Начинаем с k=3, чтобы вычислить y_4, используя историю 0,1,2,3.
        for(int k=3; k<n; ++k) {
            // Предиктор
            // --- ЭТАП ПРЕДИКТОРА ---
            // Делаем "предсказание" y_{k+1} и z_{k+1} с помощью явной формулы Адамса.
            // Это в точности формула (4.26) из методички.
            // y_p = y_k + (h/24)*(55f_k - 59f_{k-1} + ...)
            double y_p = yh.Get(k) + (h/24.0)*(55*fyh.Get(k) - 59*fyh.Get(k-1) + 37*fyh.Get(k-2) - 9*fyh.Get(k-3));
            double z_p = zh.Get(k) + (h/24.0)*(55*fzh.Get(k) - 59*fzh.Get(k-1) + 37*fzh.Get(k-2) - 9*fzh.Get(k-3));
            double next_x = xh.Get(k) + h;

            // ВАЖНО: Считаем ОБЕ производные от предсказанных, согласованных (y_p, z_p)
            // Вычисляем "предсказанные" значения производных в точке x_{k+1},
            // используя предсказанные y_p и z_p.
            // f_p = f(x_{k+1}, y_p, z_p)
            double fy_p = fy(next_x, y_p, z_p);
            double fz_p = fz(next_x, y_p, z_p);

            // Корректор
            // --- ЭТАП КОРРЕКТОРА ---
            // Делаем "уточнение" y_{k+1} и z_{k+1}, используя неявную формулу Адамса.
            // В нее подставляется "предсказанная" производная f_p.
            // Это в точности формула (4.27) из методички.
            // y_c = y_k + (h/24)*(9f_p + 19f_k - 5f_{k-1} + ...)
            double y_c = yh.Get(k) + (h/24.0)*(9*fy_p + 19*fyh.Get(k) - 5*fyh.Get(k-1) + fyh.Get(k-2));
            double z_c = zh.Get(k) + (h/24.0)*(9*fz_p + 19*fzh.Get(k) - 5*fzh.Get(k-1) + fzh.Get(k-2));

            // Сохраняем в историю
            // Сохраняем в историю СКОРРЕКТИРОВАННЫЕ (уточненные) значения.
            xh.Set(k+1, next_x);
            yh.Set(k+1, y_c);
            zh.Set(k+1, z_c);

            // Вычисляем и сохраняем производные от СКОРРЕКТИРОВАННЫХ значений.
            fyh.Set(k+1, fy(next_x, y_c, z_c));
            fzh.Set(k+1, fz(next_x, y_c, z_c));
            
            PrintTableRow(k+1, next_x, y_c, z_c);
        }
    }

    // =========================================================================
    // НОВЫЙ РАЗДЕЛ: "ЧИСТЫЕ" РЕШАТЕЛИ (только вычисления, без печати)
    // =========================================================================
    
    double SolverEulerExplicit(double step) {
        int n = static_cast<int>(round((xk - x0) / step)); double y = y0, z = z0;
        for (int k = 0; k < n; ++k) { double x = x0 + k * step; y += step * fy(x,y,z); z += step * fz(x,y,z); }
        return y;
    }

    double SolverEulerImplicit(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); double y=y0,z=z0;
        for (int k = 0; k < n; ++k) { double x=x0+k*step; double next_x=x+step; double y_next=y+step*fy(x,y,z), z_next=z+step*fz(x,y,z);
            for(int i=0;i<10;++i){ y_next=y+step*fy(next_x,y_next,z_next); z_next=z+step*fz(next_x,y_next,z_next); }
            y=y_next; z=z_next;
        } return y;
    }

    double SolverEulerCauchy(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); double y=y0,z=z0;
        for (int k = 0; k < n; ++k) { double x=x0+k*step; double next_x=x+step; double y_pred=y+step*fy(x,y,z), z_pred=z+step*fz(x,y,z);
            y+=(step/2.0)*(fy(x,y,z)+fy(next_x,y_pred,z_pred)); z+=(step/2.0)*(fz(x,y,z)+fz(next_x,y_pred,z_pred));
        } return y;
    }

    double SolverEulerCauchyIterative(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); double y=y0,z=z0;
        for (int k = 0; k < n; ++k) { double x=x0+k*step; double next_x=x+step; double y_next=y+step*fy(x,y,z), z_next=z+step*fz(x,y,z);
            double fy0=fy(x,y,z), fz0=fz(x,y,z);
            for(int i=0; i<10; ++i){ y_next=y+(step/2.0)*(fy0+fy(next_x,y_next,z_next)); z_next=z+(step/2.0)*(fz0+fz(next_x,y_next,z_next)); }
            y=y_next; z=z_next;
        } return y;
    }

    double SolverEulerImprovedMidpoint(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); double y=y0,z=z0;
        for (int k = 0; k < n; ++k) { double x=x0+k*step; double mid_x=x+step/2.0; double y_mid=y+(step/2.0)*fy(x,y,z), z_mid=z+(step/2.0)*fz(x,y,z);
            y+=step*fy(mid_x,y_mid,z_mid); z+=step*fz(mid_x,y_mid,z_mid);
        } return y;
    }

    double SolverRK3(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); double y=y0,z=z0;
        for (int k = 0; k < n; ++k) { double x=x0+k*step;
            double K1=step*fy(x,y,z), L1=step*fz(x,y,z); double K2=step*fy(x+step/3.,y+K1/3.,z+L1/3.), L2=step*fz(x+step/3.,y+K1/3.,z+L1/3.);
            double K3=step*fy(x+2.*step/3.,y+2.*K2/3.,z+2.*L2/3.), L3=step*fz(x+2.*step/3.,y+2.*K2/3.,z+2.*L2/3.);
            y+=(K1+3.*K3)/4.; z+=(L1+3.*L3)/4.;
        } return y;
    }

    double SolverRK4(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); double y=y0,z=z0;
        for (int k=0; k<n; ++k){ double x=x0+k*step;
            double K1=step*fy(x,y,z), L1=step*fz(x,y,z), K2=step*fy(x+step/2,y+K1/2,z+L1/2), L2=step*fz(x+step/2,y+K1/2,z+L1/2);
            double K3=step*fy(x+step/2,y+K2/2,z+L2/2), L3=step*fz(x+step/2,y+K2/2,z+L2/2); double K4=step*fy(x+step,y+K3,z+L3), L4=step*fz(x+step,y+K3,z+L3);
            y+=(K1+2*K2+2*K3+K4)/6.0; z+=(L1+2*L2+2*L3+L4)/6.0;
        } return y;
    }
    
    double SolverAdamsExplicit(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); if(n<4) return NAN;
        Vector xh(n+1), yh(n+1), zh(n+1), fyh(n+1), fzh(n+1); RungeKuttaStartup(xh,yh,zh,fyh,fzh,step);
        for (int k=3; k<n; ++k) {
            double y_n=yh.Get(k)+(step/24.)*(55*fyh.Get(k)-59*fyh.Get(k-1)+37*fyh.Get(k-2)-9*fyh.Get(k-3));
            double z_n=zh.Get(k)+(step/24.)*(55*fzh.Get(k)-59*fzh.Get(k-1)+37*fzh.Get(k-2)-9*fzh.Get(k-3));
            double x_n=xh.Get(k-1)+step; yh.Set(k+1,y_n); zh.Set(k+1,z_n); xh.Set(k+1,x_n);
            fyh.Set(k+1,fy(x_n,y_n,z_n)); fzh.Set(k+1,fz(x_n,y_n,z_n));
        } return yh.Get(n);
    }
    
    double SolverAdamsPredictorCorrector(double step) {
        int n=static_cast<int>(round((xk-x0)/step)); if(n<4) return NAN;
        Vector xh(n+1), yh(n+1), zh(n+1), fyh(n+1), fzh(n+1); RungeKuttaStartup(xh,yh,zh,fyh,fzh,step);
        for (int k=3; k<n; ++k) {
            double y_p=yh.Get(k)+(step/24.)*(55*fyh.Get(k)-59*fyh.Get(k-1)+37*fyh.Get(k-2)-9*fyh.Get(k-3));
            double z_p=zh.Get(k)+(step/24.)*(55*fzh.Get(k)-59*fzh.Get(k-1)+37*fzh.Get(k-2)-9*fzh.Get(k-3));
            double x_n=xh.Get(k)+step; double fy_p=fy(x_n,y_p,z_p), fz_p=fz(x_n,y_p,z_p);
            double y_c=yh.Get(k)+(step/24.)*(9*fy_p+19*fyh.Get(k)-5*fyh.Get(k-1)+fyh.Get(k-2));
            double z_c=zh.Get(k)+(step/24.)*(9*fz_p+19*fzh.Get(k)-5*fzh.Get(k-1)+fzh.Get(k-2));
            xh.Set(k+1,x_n); yh.Set(k+1,y_c); zh.Set(k+1,z_c);
            fyh.Set(k+1,fy(x_n,y_c,z_c)); fzh.Set(k+1,fz(x_n,y_c,z_c));
        } return yh.Get(n);
    }

    // =========================================================================
    // ОБНОВЛЕННЫЙ БЛОК ОЦЕНКИ ПОГРЕШНОСТИ
    // =========================================================================

    // НОВАЯ универсальная функция для вычисления и печати ошибок для ЛЮБОГО метода
    void PrintFormattedErrorLine(const std::string& name, const std::string& val1, const std::string& val2) {
        // Задаем ширину колонок в символах
        const int name_width = 38;
        const int error_width = 25;

        std::string line = name;
        // Дополняем имя метода пробелами справа до нужной ширины
        if (name.length() < name_width) {
            line += std::string(name_width - name.length(), ' ');
        }
        
        std::string val1_padded = val1;
        // Дополняем первое число пробелами слева до нужной ширины (выравнивание по правому краю)
        if (val1.length() < error_width) {
            val1_padded = std::string(error_width - val1.length(), ' ') + val1;
        }

        std::string val2_padded = val2;
        // Дополняем второе число пробелами слева до нужной ширины
        if (val2.length() < error_width) {
            val2_padded = std::string(error_width - val2.length(), ' ') + val2;
        }
        
        std::cout << line << val1_padded << val2_padded << "\n";
    }

    // ОСНОВНАЯ функция для оценки погрешности.
    void EvaluateRungeRombergForAll() {
        std::cout << "\n\n--- Оценка погрешности в конечной точке x=" << xk << " ---\n";

        // Создаем лямбда-функцию для вывода одного блока, чтобы не повторять код
        auto PrintResultBlock = [&](const std::string& name, int p, std::function<double(double)> solver) {

            // ---= ШАГ 2: Двойное вычисление (суть метода Рунге-Ромберга) =---
            // Вызываем "тихий" решатель solver (например, SolverRK4) с шагом h (например, 0.1)
            double y_h = solver(h);

            // Вызываем тот же решатель, но с удвоенным шагом 2*h (например, 0.2)
            double y_2h = solver(h * 2.0);

            std::cout << "\n=======================================================\n";
            std::cout << "[ " << name << " ]\n";
            std::cout << "-------------------------------------------------------\n";

            // Используем stringstream для преобразования чисел в строки с нужной точностью
            std::stringstream ss_rr, ss_true;

            if (std::isnan(y_h) || std::isnan(y_2h)) {
                std::cout << "  Недостаточно точек для оценки.\n";
            } else {

                // ---= ШАГ 3: Вычисление ДВУХ видов погрешности =---

                // 3.1. Оценка по Рунге-Ромбергу (Практическая).
                //      Это в точности формула (4.11) из методички: (y_h - y_2h) / (2^p - 1)
                //      p - порядок точности метода (1 для Эйлера, 4 для РК4 и Адамса).
                double error_rr = std::abs(y_h - y_2h) / (pow(2, p) - 1.0);

                // 3.2. Сравнение с точным решением (Истинная ошибка).
                //      Это просто разница между нашим лучшим численным результатом (y_h)
                //      и значением из точной аналитической формулы.
                double error_true = std::abs(y_exact(xk) - y_h);
                
                ss_rr << std::fixed << std::setprecision(8) << error_rr;
                ss_true << std::fixed << std::setprecision(8) << error_true;

                // Выводим с выравниванием только для двоеточий
                std::cout << std::left << std::setw(25) << "  Оценка по Р-Р" << ": " << ss_rr.str() << "\n";
                std::cout << std::left << std::setw(25) << "  Истинная ошибка" << ": " << ss_true.str() << "\n";
            }
            std::cout << "=======================================================\n";
        };

        // Последовательно вызываем для каждого метода
        // ---= ШАГ 4: Последовательный вызов "шаблона" для каждого метода =---
        // Теперь мы просто 9 раз вызываем наш PrintResultBlock, каждый раз передавая ему:
        // 1. Имя метода для печати (например, "1. Явный метод Эйлера")
        // 2. Порядок точности 'p' для этого метода (например, 1)
        // 3. Сам "тихий" решатель, который нужно использовать (например, SolverEulerExplicit)
        PrintResultBlock("1. Явный метод Эйлера", 1, [&](double s){ return SolverEulerExplicit(s); });
        PrintResultBlock("2. Неявный метод Эйлера", 1, [&](double s){ return SolverEulerImplicit(s); });
        PrintResultBlock("3. Эйлер-Коши (явный)", 2, [&](double s){ return SolverEulerCauchy(s); });
        PrintResultBlock("4. Эйлер-Коши (неявный)", 2, [&](double s){ return SolverEulerCauchyIterative(s); });
        PrintResultBlock("5. Улучшенный Эйлер (ср. точка)", 2, [&](double s){ return SolverEulerImprovedMidpoint(s); });
        PrintResultBlock("6. Рунге-Кутта 3-го порядка", 3, [&](double s){ return SolverRK3(s); });
        PrintResultBlock("7. Рунге-Кутта 4-го порядка", 4, [&](double s){ return SolverRK4(s); });
        PrintResultBlock("8. Адамс явный 4-го порядка", 4, [&](double s){ return SolverAdamsExplicit(s); });
        PrintResultBlock("9. Адамс П-К 4-го порядка", 4, [&](double s){ return SolverAdamsPredictorCorrector(s); });
    }
    
    // Вспомогательные функции (без изменений)
    void PrintTableHeader(const std::string& title) {
        std::cout << "\n\n--- " << title << " ---\n";
        std::cout << std::setw(5) << "k" << std::setw(12) << "x_k" << std::setw(15) << "y_k" << std::setw(15) << "z_k"
                  << std::setw(15) << "dy_k" << std::setw(15) << "dz_k" // <-- ДОБАВЛЕНЫ СТОЛБЦЫ
                  << std::setw(15) << "y_ист" << std::setw(15) << "epsilon_k" << "\n";
        std::cout << "------------------------------------------------------------------------------------------------------\n"; // Линия стала длиннее
    }

    void PrintTableRow(int k, double x, double y, double z) {
        double exact_y = y_exact(x);
        double error = std::abs(exact_y - y);
        std::cout << std::setw(5) << k << std::setw(12) << x << std::setw(15) << y << std::setw(15) << z
                  << std::setw(15) << " " << std::setw(15) << " " // Пустые поля для dy, dz
                  << std::setw(15) << exact_y << std::setw(15) << error << "\n";
    }

    void PrintTableRow(int k, double x, double y, double z, double dy, double dz) {
        double exact_y = y_exact(x);
        double error = std::abs(exact_y - y);
        std::cout << std::setw(5) << k << std::setw(12) << x << std::setw(15) << y << std::setw(15) << z
                  << std::setw(15) << dy << std::setw(15) << dz // <-- ВЫВОДИМ dy и dz
                  << std::setw(15) << exact_y << std::setw(15) << error << "\n";
    }

    // вычисляем первые 4 точки с помощью Рунге-Кутты 4-го порядка
    void RungeKuttaStartup(Vector& xh, Vector& yh, Vector& zh, Vector& fyh, Vector& fzh, double step) {
        double cx=x0, cy=y0, cz=z0;

        // 2. Делаем 4 итерации, чтобы заполнить историю для k=0, 1, 2, 3
        for (int k=0; k<4; ++k) {

            // 3. Сохраняем ТЕКУЩИЕ значения в "историю".
            //    На первой итерации (k=0) сюда запишутся начальные x0, y0, z0.
            xh.Set(k, cx); yh.Set(k, cy); zh.Set(k, cz);

            // 4. Сразу вычисляем и сохраняем производные в текущей точке.
            fyh.Set(k, fy(cx,cy,cz)); fzh.Set(k, fz(cx,cy,cz));

            // 5. Если мы уже заполнили последнюю, 4-ю точку (k=3), выходим из цикла.
            //    Дальше считать не нужно, "разгон" завершен.
            if (k==3) break;
            
            // 6. Вычисляем ОДИН шаг методом Рунге-Кутты 4-го порядка.

            // ---= РЕАЛИЗАЦИЯ ОДНОГО ШАГА МЕТОДА РУНГЕ-КУТТЫ 4-го ПОРЯДКА (Формулы 4.10, 4.16) =---
            // K_i - коэффициенты для уравнения y' = fy(...)
            // L_i - коэффициенты для уравнения z' = fz(...)

            // Шаг 1: Вычисляем коэффициенты K1, L1 на основе наклона в НАЧАЛЬНОЙ точке.
            // Формула: K1 = h * f(x_k, y_k, z_k)

            // Шаг 2: Вычисляем коэффициенты K2, L2 на основе наклона в СРЕДНЕЙ точке,
            //        которую мы предсказываем с помощью K1, L1.
            // Формула: K2 = h * f(x_k + h/2, y_k + K1/2, z_k + L1/2)

            // Шаг 3: Вычисляем коэффициенты K3, L3 на основе УТОЧНЕННОГО наклона в СРЕДНЕЙ точке,
            //        используя более точные K2, L2.
            // Формула: K3 = h * f(x_k + h/2, y_k + K2/2, z_k + L2/2)

            // Шаг 4: Вычисляем коэффициенты K4, L4 на основе наклона в КОНЕЧНОЙ точке,
            //        которую мы предсказываем с помощью K3, L3.
            // Формула: K4 = h * f(x_k + h, y_k + K3, z_k + L3)
            double K1=step*fy(cx,cy,cz), L1=step*fz(cx,cy,cz), K2=step*fy(cx+step/2,cy+K1/2,cz+L1/2), L2=step*fz(cx+step/2,cy+K1/2,cz+L1/2);
            double K3=step*fy(cx+step/2,cy+K2/2,cz+L2/2), L3=step*fz(cx+step/2,cy+K2/2,cz+L2/2), K4=step*fy(cx+step,cy+K3,cz+L3), L4=step*fz(cx+step,cy+K3,cz+L3);

            // 7. Обновляем "текущие" переменные (cx, cy, cz) до значений
        //    на следующем шаге, чтобы использовать их в следующей итерации цикла.
            cy += (K1+2*K2+2*K3+K4)/6.0; cz += (L1+2*L2+2*L3+L4)/6.0; cx += step;
        }
    }
    
    // Эту функцию нужно будет менять в main в зависимости от задачи
    void PrintInitialData() const {
        std::cout << "--- Исходные данные (Вариант 13) ---\n";
        std::cout << "ОДУ: y'' = 2*tg(x)*y' + 3*y\n";
        std::cout << "Начальные условия: y(0) = 1, y'(0) = 3\n";
        std::cout << "Отрезок: [" << x0 << ", " << xk << "], Шаг h = " << h << "\n";
        std::cout << "Точное решение: y = cos^3(x) + sin(x)*(1 + 2*cos^2(x))\n";
    }
};

// 1 вариант
// int main() {
//     try {
//         // y' = z
//         ODEFunc func_y = [](double, double, double z) { return z; };
//         // y'' = sin(3x) - y
//         ODEFunc func_z = [](double x, double y, double) { return sin(3.0 * x) - y; };
        
//         // y = cos(x) + (11/8)*sin(x) - sin(3x)/8
//         ExactSolutionFunc exact_sol = [](double x) { 
//             return cos(x) + (11.0 / 8.0) * sin(x) - sin(3.0 * x) / 8.0; 
//         };

//         // Начальные условия: y(0)=1, y'(0)=1
//         double x_start = 0.0, y_start = 1.0, z_start = 1.0, x_end = 1.0, h_step = 0.1;

//         Task_4_1 task(func_y, func_z, exact_sol, x_start, y_start, z_start, x_end, h_step);
        
//         // Не забудьте обновить текстовое описание, если хотите
//         // Это можно сделать в методе PrintInitialData()
//         task.Do();

//     } catch (const std::exception& e) {
//         std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
//         return 1;
//     }
//     return 0;
// }

int main() {
    try {
        ODEFunc func_y = [](double, double, double z) { return z; };
        // ODEFunc func_z = [](double x, double y, double z) { return 2.0 * tan(x) * z + 3.0 * y; }; //!!!!
        ODEFunc func_z = [](double x, double y, double z) { return -2.0 * tan(x) * z - 3.0 * y; };
        ExactSolutionFunc exact_sol = [](double x) { return pow(cos(x), 3) + sin(x) * (1.0 + 2.0 * pow(cos(x), 2)); };

        double x_start = 0.0, y_start = 1.0, z_start = 3.0, x_end = 1.0, h_step = 0.1;

        Task_4_1 task(func_y, func_z, exact_sol, x_start, y_start, z_start, x_end, h_step);
        task.Do();
    } catch (const std::exception& e) {
        std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

//Методичка 1 функция
// int main() {
//     try {
//         // Наша основная функция: y' = (y+x)^2
//         // Параметр z игнорируется, так как уравнение первого порядка.
//         ODEFunc func_y = [](double x, double y, double /*z_ignored*/) {
//             return (y + x) * (y + x);
//         };
        
//         // Так как уравнение первого порядка, второй функции z у нас нет.
//         // Создаем "пустышку", которая всегда возвращает 0.
//         ODEFunc func_z = [](double, double, double) {
//             return 0.0;
//         };
        
//         ExactSolutionFunc exact_sol = [](double x) { 
//             return tan(x) - x; 
//         };

//         double x_start = 0.0;
//         double y_start = 0.0;
//         double z_start = 0.0; // Начальное условие для нашей фиктивной переменной z
//         // Интервал [0, 0.5] с шагом h=0.1
//         double x_end = 0.5;
//         double h_step = 0.1;

//         // 4. Создаем и запускаем задачу
//         Task_4_1 task(func_y, func_z, exact_sol, x_start, y_start, z_start, x_end, h_step);
        
//         // ВАЖНО: нужно будет обновить метод PrintInitialData() внутри класса,
//         // чтобы он выводил правильное описание для этой новой задачи.
//         task.Do();

//     } catch (const std::exception& e) {
//         std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
//         return 1;
//     }
//     return 0;
// }

// Методичка 2 функция
// int main() {
//     try {
//         // --- Настройка для Примера 4.5 из методички ---

//         // 1. Сводим уравнение (x^2+1)y'' = 2xy' к системе.
//         //    y'' = (2*x*y') / (x^2 + 1)
//         //    Вводим замену z = y', получаем систему:
//         //    y' = z
//         //    z' = (2*x*z) / (x^2 + 1)

//         // 2. Определяем правые части системы
//         ODEFunc func_y = [](double x, double y, double z) {
//             return z; // y' = z
//         };
//         ODEFunc func_z = [](double x, double y, double z) {
//             return (2.0 * x * z) / (x * x + 1.0); // z' = (2*x*z)/(x^2+1)
//         };
        
//         // 3. Точное решение из примера
//         ExactSolutionFunc exact_sol = [](double x) { 
//             return pow(x, 3) + 3.0 * x + 1.0; 
//         };

//         // 4. Параметры задачи из примера
//         // Начальные условия: y(0) = 1, y'(0) = 3
//         double x_start = 0.0;
//         double y_start = 1.0;
//         double z_start = 3.0;
//         // Интервал [0, 1] с шагом h=0.2
//         double x_end = 1.0;
//         double h_step = 0.2;

//         // 5. Создаем и запускаем задачу
//         Task_4_1 task(func_y, func_z, exact_sol, x_start, y_start, z_start, x_end, h_step);
        
//         // ВАЖНО: Не забудь обновить метод PrintInitialData() внутри класса,
//         // чтобы он выводил правильное описание для этой задачи.
//         task.Do();

//     } catch (const std::exception& e) {
//         std::cerr << "\nПроизошла ошибка: " << e.what() << std::endl;
//         return 1;
//     }
//     return 0;
// }