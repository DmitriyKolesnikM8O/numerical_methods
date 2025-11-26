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
            double dy = h * fy(x,y,z);
            double dz = h * fz(x,y,z);
            PrintTableRow(k, x, y, z, dy, dz);
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
            double y_next = y + h*fy(x,y,z);
            double z_next = z + h*fz(x,y,z);
            // Итерации для поиска y_{k+1} и z_{k+1}
            for(int i=0;i<10;++i){
                double y_iter_old = y_next;
                y_next = y + h * fy(next_x, y_next, z_next);
                z_next = z + h * fz(next_x, y_iter_old, z_next);
            }

            // ВЫЧИСЛЯЕМ ПРИРАЩЕНИЯ ПОСТФАКТУМ
            double dy = y_next - y;
            double dz = z_next - z;
            PrintTableRow(k, x, y, z, dy, dz); // Печатаем с приращениями
            
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
            double y_pred = y + h * fy(x,y,z), z_pred = z + h * fz(x,y,z);
            double dy = (h/2.0) * (fy(x,y,z) + fy(next_x, y_pred, z_pred));
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
            double y_next = y + h*fy(x,y,z);
            double z_next = z + h*fz(x,y,z);
            double fy0 = fy(x,y,z), fz0 = fz(x,y,z);
            // Итерации для поиска y_{k+1} и z_{k+1}
            for(int i=0; i<10; ++i){
                double y_iter_old = y_next;
                y_next = y + (h/2.0) * (fy0 + fy(next_x, y_next, z_next));
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
            double mid_x = x + h/2.0;
            double y_mid = y + (h/2.0) * fy(x,y,z);
            double z_mid = z + (h/2.0) * fz(x,y,z);
            double dy = h * fy(mid_x, y_mid, z_mid);
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
            double K1=h*fy(x,y,z), L1=h*fz(x,y,z);
            double K2=h*fy(x+h/3.,y+K1/3.,z+L1/3.), L2=h*fz(x+h/3.,y+K1/3.,z+L1/3.);
            double K3=h*fy(x+2.*h/3.,y+2.*K2/3.,z+2.*L2/3.), L3=h*fz(x+2.*h/3.,y+2.*K2/3.,z+2.*L2/3.);
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
            double K1=h*fy(x,y,z), L1=h*fz(x,y,z);
            double K2=h*fy(x+h/2,y+K1/2,z+L1/2), L2=h*fz(x+h/2,y+K1/2,z+L1/2);
            double K3=h*fy(x+h/2,y+K2/2,z+L2/2), L3=h*fz(x+h/2,y+K2/2,z+L2/2);
            double K4=h*fy(x+h,y+K3,z+L3), L4=h*fz(x+h,y+K3,z+L3);
            
            // Вычисляем ПРИРАЩЕНИЯ dy и dz
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
        PrintTableHeader("9. Метод Адамса (явный, 4-го порядка)"); // Используем PrintTableRow без dy, dz
        int n=static_cast<int>(round((xk-x0)/h)); if(n<4){std::cout<<"Мало шагов\n"; return;}
        Vector xh(n+1), yh(n+1), zh(n+1), fyh(n+1), fzh(n+1);
        RungeKuttaStartup(xh, yh, zh, fyh, fzh, h);

        for(int k=0; k<4; ++k) PrintTableRow(k, xh.Get(k), yh.Get(k), zh.Get(k));

        for(int k=3; k<n; ++k) {
            // Получаем значения y_k и z_k из истории
            double current_y = yh.Get(k);
            double current_z = zh.Get(k);

            // Считаем y_{k+1} и z_{k+1}
            double next_y = current_y + (h/24.0) * (55*fyh.Get(k) - 59*fyh.Get(k-1) + 37*fyh.Get(k-2) - 9*fyh.Get(k-3));
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
        PrintTableHeader("10. Метод Адамса-Бэшфортса-Моултона");
        int n=static_cast<int>(round((xk-x0)/h)); if(n<4){std::cout<<"Мало шагов\n"; return;}
        Vector xh(n+1), yh(n+1), zh(n+1), fyh(n+1), fzh(n+1);
        RungeKuttaStartup(xh, yh, zh, fyh, fzh, h);

        for(int k=0; k<4; ++k) PrintTableRow(k, xh.Get(k), yh.Get(k), zh.Get(k));
        
        for(int k=3; k<n; ++k) {
            // Предиктор
            double y_p = yh.Get(k) + (h/24.0)*(55*fyh.Get(k) - 59*fyh.Get(k-1) + 37*fyh.Get(k-2) - 9*fyh.Get(k-3));
            double z_p = zh.Get(k) + (h/24.0)*(55*fzh.Get(k) - 59*fzh.Get(k-1) + 37*fzh.Get(k-2) - 9*fzh.Get(k-3));
            double next_x = xh.Get(k) + h;

            // ВАЖНО: Считаем ОБЕ производные от предсказанных, согласованных (y_p, z_p)
            double fy_p = fy(next_x, y_p, z_p);
            double fz_p = fz(next_x, y_p, z_p);

            // Корректор
            double y_c = yh.Get(k) + (h/24.0)*(9*fy_p + 19*fyh.Get(k) - 5*fyh.Get(k-1) + fyh.Get(k-2));
            double z_c = zh.Get(k) + (h/24.0)*(9*fz_p + 19*fzh.Get(k) - 5*fzh.Get(k-1) + fzh.Get(k-2));

            // Сохраняем в историю
            xh.Set(k+1, next_x);
            yh.Set(k+1, y_c);
            zh.Set(k+1, z_c);
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
            
            double y_h = solver(h);
            double y_2h = solver(h * 2.0);

            std::cout << "\n=======================================================\n";
            std::cout << "[ " << name << " ]\n";
            std::cout << "-------------------------------------------------------\n";

            // Используем stringstream для преобразования чисел в строки с нужной точностью
            std::stringstream ss_rr, ss_true;

            if (std::isnan(y_h) || std::isnan(y_2h)) {
                std::cout << "  Недостаточно точек для оценки.\n";
            } else {
                double error_rr = std::abs(y_h - y_2h) / (pow(2, p) - 1.0);
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

    void RungeKuttaStartup(Vector& xh, Vector& yh, Vector& zh, Vector& fyh, Vector& fzh, double step) {
        double cx=x0, cy=y0, cz=z0;
        for (int k=0; k<4; ++k) {
            xh.Set(k, cx); yh.Set(k, cy); zh.Set(k, cz);
            fyh.Set(k, fy(cx,cy,cz)); fzh.Set(k, fz(cx,cy,cz));
            if (k==3) break;
            double K1=step*fy(cx,cy,cz), L1=step*fz(cx,cy,cz), K2=step*fy(cx+step/2,cy+K1/2,cz+L1/2), L2=step*fz(cx+step/2,cy+K1/2,cz+L1/2);
            double K3=step*fy(cx+step/2,cy+K2/2,cz+L2/2), L3=step*fz(cx+step/2,cy+K2/2,cz+L2/2), K4=step*fy(cx+step,cy+K3,cz+L3), L4=step*fz(cx+step,cy+K3,cz+L3);
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