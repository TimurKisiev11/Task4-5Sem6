// Task4Sem6.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <windows.h>
#include <iomanip>
#include <cstdlib>
using namespace std;
//выводит матриицу A порядка n
void show_mat(double** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << fixed << setprecision(6) << " \t" << A[i][j] << " \t";
        }
        std::cout << endl;
    }
}
//строит нулевой фектор порядка n
double* get_null_vect(int n)
{
    double* V = new double[2 * n];
    for (int i = 0; i <= n; i++)
    {
        V[i] = 0.0;
    }
    return V;
}
//строит нулевую матрицу порядка n
double** get_null_mat(int n)
{
    double** A = new double* [2 * n];
    for (int i = 0; i <= n; i++)
    {
        A[i] = new double[2 * n];
        for (int j = 0; j <= n; j++)
            A[i][j] = 0.0;
    }
    return A;
}
//копирует матрицу B в A
void cop(double** a, double** b, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i][j] = b[i][j];
}
//строит нулевую матрицу 
double** get_null_matrix(int n)
{
    n = n + 1;
    double** A = new double* [2 * n];
    for (int i = 0; i <= n; i++)
    {
        A[i] = new double[2 * n];
        for (int j = 0; j <= n; j++)
            A[i][j] = 0.0;
    }
    return A;
}

//Многочлены Лежандра (для нахождения узлов для формулы Гауса)
double Legendre(int n, double x)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return x;
    else
        return (Legendre(n - 1, x) * x * (2 * n - 1) / (n)) - (Legendre(n - 2, x) * (n - 1) / (n));
}
//Локализация корней многочленов Лежандра  (для нахождения узлов для формулы Гауса)
double* loc_root(double a, double b, double N, int n)
{
    double* mas = new double[100];
    mas[0] = n;
    double h = (b - a) / N;
    double tmp = 0;
    int j = 2;
    while (a + h <= b)
    {
        if ((Legendre(n, a) * Legendre(n, a + h)) <= 0)
        {
            mas[j] = a;
            mas[j + 1] = a + h;
            j += 2;
        }
        a = a + h;
    }
    if (mas[0] > 0)
    {
        //cout << fixed << setprecision(0) << "количество корней ортоганального многочлена на отрезке (a,b) = " << mas[0] << endl;
        return mas;
    }
    else
    {
        //cout << "нет корней на данном отрезке" << endl;
    }
    return 0;
}

//Коэффициенты для формулы Гаусса
double* koeff(double* z, int n)
{
    double* A = new double[100];
    for (int i = 0; i < n; i++)
    {
        A[i] = (2 * (1 - pow(z[i], 2)) / (pow(n, 2) * Legendre(n - 1, z[i]) * Legendre(n - 1, z[i])));
    }
    return A;
}

//Поиск корней многочленов Лежандра  (для нахождения узлов для формулы Гауса)
double Secant(double a, double b, double EPS, int n)
{
    double fapprox, sapprox, x;
    int counter = 0;

    fapprox = a;
    sapprox = b;

    x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
    counter++;

    while (abs(x - sapprox) >= EPS)
    {
        fapprox = sapprox;
        sapprox = x;
        x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
        counter++;
    }

    return x;
}

//Вернет матрицу для каждого отрезка [a,b]; Первая строка - узлы КФ Гаусса, вторая - коэффициенты.
double** Gauss(double a, double b, int n)
{
    double** result = get_null_mat(2 * n);
    double* z = get_null_vect(2 * n);        //массив для узлов для формулы Гаусса
    int val = n;                             // количество узлов для формулы Гаусса
    double* MASS = loc_root(-1, 1, 100, val);//локализуем узлы многочлена Лежандра
    double counter = MASS[0];
    int cnt = 0;
    for (int j = 2; j <= 2 * counter + 1; j += 2)
    {
        z[cnt] = Secant(MASS[j], MASS[j + 1], 0.0000000001, val);//находим узлы (корни многочлена Лежандра)
        cnt++;
    }
    if (val % 2 > 0)
    {
        int m = val / 2;
        z[m] = 0.0;
    }
    double* A = koeff(z, val); //находим коэффициенты формулы гаусса
    // cout << "Узлы и коэффициенты для формулы Гаусса: " << endl;
    // cout << "      Узел      " << " <-> " << "    Коэффициент    " << endl;
    for (int k = 0; k < cnt; k++)
    {
        z[k] = (b - a) / 2.0 * z[k] + (b + a) / 2.0; //узлы для формулы гаусса
        A[k] = (b - a) / 2.0 * A[k];
       // cout << fixed << setprecision(2) << z[k] << "            <->            " << A[k] << endl;
    }
    result[0] = z;
    result[1] = A;
    // cout << "Соберем в матрицу узлы и коэффициенты" << endl;
    // show_mat(result, n);
    return result;
}
//функция phi из условия 
double phi(double x)
{
    double a = 2 * M_PI * x;
    double b = 3 * M_PI * x;
    return(2 * sin(a) + (-1) * sin(b));
}
//ортонормированные собственные функции
double psi(double p, double x)
{
    double a = p * M_PI * x;
    return(sqrt(2) * sin(a));
}
//получает точное решение в точке (x,t)
double accurate(int N, double x, double t)
{
    double** gauss = Gauss(0, 1, 9);
    double* c_p = new double[2 * N];
    double u_x_t = 0.0;
    for (int p = 1; p < N; p++) 
    {
        c_p[p] = 0.0;
        for (int k = 0; k < 9; k++)
        {
            c_p[p] += gauss[1][k] * (phi(gauss[0][k])*psi(p,gauss[0][k]));
        }
       // cout << p << " -- " << c_p[p] << endl;
    }
    for (int p = 1; p < N; p++)
    {
        double exp_arg = (-1) * pow(M_PI, 2) * pow(p, 2) * t;
        double sin_arg = p * M_PI * x;
        u_x_t += (c_p[p] * exp(exp_arg) * psi(p,x));
       // cout<<fixed<<setprecision(20) << "  " << p << " " << u_x_t << endl;
    }
    return ( u_x_t);
}
//находит решение в точке (x,t) ДРФ
double DFT(int N, double x, double t)
{
    double h = 1.0 / N;
    double* c_p = new double[2 * N];
    double u_x_t = 0.0;
    for (int p = 1; p < N; p++)
    {
        c_p[p] = 0.0;
        for (int i = 1; i < N; i++)
        {
            //cout << phi(i * h) << "   ----    " << sin(p * M_PI * i * h) << endl;
            c_p[p] += phi(i * h) * sin(p * M_PI * i * h);
        }
        c_p[p] *= (sqrt(2) * h);
       // cout << "  " << p << " " << c_p[p] << endl;
    }
    for (int p = 1; p < N; p++)
    {   
        double exp_arg = (-1) * pow(M_PI, 2) * pow(p, 2) * t;
        double sin_arg = p * M_PI * x;
        u_x_t += (c_p[p] * exp(exp_arg) * psi(p, x));
        //cout << "  " << p << " " << u_x_t << endl;
    }
   // cout <<sqrt(2)* u_x_t << endl;
    return (u_x_t);
}
/*находит решение в точке(x, t) ДРФ с весами
flag = 1 => sigma = 0
flag = 2 => sigma = 1
flag = 3 => sigma = 0.5
flag = 4 => sigma = 0.5 - (pow(h, 2) / (12 * tau))
*/
double DFT_NET(int flag, int N, double x, double t)
{
    double h = 1.0 / N;
    double T = 0.1;
    int M =  10 * N;
    double tau = T / M;
    double sigma = 0.0;
    double* lambda = new double[2 * N];
    double r = tau / pow(h, 2);
    double sin_arg = 0.0;
    switch (flag)
    {
    case 1:
    { sigma = 0.0; }
    break;
    case 2:
    {sigma = 1.0; }
    break;
    case 3:
    { sigma = 0.5; }
    break;
    case 4:
    {sigma = 0.5 - (pow(h, 2) / (12 * tau)); }
    break;
    default:
    {
        sigma = 0.0;
    }
    break;
    }
    for (int p = 1; p < N; p++)
    {   
       
        sin_arg = p * M_PI * h;
        double numinator = 1 - (4 * (1 - sigma) * r * pow(sin(sin_arg / 2.0), 2));
        double denominator = 1 + (4 * sigma * r * pow(sin(sin_arg / 2.0), 2));
        lambda[p] = numinator / denominator;
       // cout << h <<" -- "<<tau<<" -- "<< pow(h,2)/(2*(1-2*sigma)) << " -- " << lambda[p] << endl;
    }
    double* c_p = new double[2 * N];
    double u_x_t = 0.0;
    for (int p = 1; p < N; p++)
    {
        c_p[p] = 0.0;
        for (int i = 1; i < N; i++)
        {
            //cout << phi(i * h) << "   ----    " << sin(p * M_PI * i * h) << endl;
            c_p[p] += h * phi(i * h) * psi(p, i * h);
        }
        // cout << "  " << p << " " << c_p[p] << endl;
    }
    double k = t / tau;
    //cout << "kkkkkkkkkkkkkk      -     " << k << endl;
    for (int p = 1; p < N; p++)
    {
        double sin_arg = p * M_PI * x;
        u_x_t += (c_p[p] * pow(lambda[p],k) * psi(p,x));
        //cout << "  " << p << " " << u_x_t << endl;
    }

    return (u_x_t);
}
//Та же функция, что выше, но выводит больше подробностей для определения устойчивости
double DFT_NET_2(int flag, int N, double x, double t)
{
    double h = 1.0 / N;
    double T = 0.1;
    int M = 5 * N;
    double tau = T / M;
    double sigma = 0.0;
    double* lambda = new double[2 * N];
    double r = tau / pow(h, 2);
    double sin_arg = 0.0;
    switch (flag)
    {
    case 1:
    { sigma = 0.0; }
    break;
    case 2:
    {sigma = 1.0; }
        break;
    case 3:
    { sigma = 0.5; }
    break;
    case 4:
    {sigma = 0.5 - (pow(h, 2) / (12 * tau)); }
    break;
    default:
    {
        sigma = 0.0;
    }
    break;
    }
    cout << " N = " << N << endl;
    cout << " sigma = " << sigma << endl;
    cout << " h = " << h << endl;
    cout << " tau = " << tau << endl;
    cout << " Собственные числа оператора перехода на следующий слой:" << endl;
    for (int p = 1; p < N; p++)
    {
        sin_arg = p * M_PI * h;
        double numinator = 1 - (4 * (1 - sigma) * r * pow(sin(sin_arg / 2.0), 2));
        double denominator = 1 + (4 * sigma * r * pow(sin(sin_arg / 2.0), 2));
        lambda[p] = numinator / denominator;
        cout <<" lambda["<<p<<"] = " << lambda[p] << endl;
    }
    if (tau <= (pow(h, 2) / (2 - 4 * sigma))) 
    {
       //cout << (pow(h, 2) / (2 - 4 * sigma)) << endl;
        cout << " Условие устойчивости для разностной схемы выполнено"<<endl;
    }
    else
    {
       //cout << (pow(h, 2) / (2 - 4 * sigma)) << endl;
        cout << " Условие устойчивости для разностной схемы НЕ выполнено" << endl;
    }
    double* c_p = new double[2 * N];
    double u_x_t = 0.0;
    for (int p = 1; p < N; p++)
    {
        c_p[p] = 0.0;
        for (int i = 1; i < N; i++)
        {
            //cout << phi(i * h) << "   ----    " << sin(p * M_PI * i * h) << endl;
            c_p[p] += h * phi(i * h) * psi(p, i * h);
        }
        // cout << "  " << p << " " << c_p[p] << endl;
    }
    double k = t / tau;
    for (int p = 1; p < N; p++)
    {
        double sin_arg = p * M_PI * x;
        u_x_t += (c_p[p] * pow(lambda[p], k) * psi(p, x));
        //cout << "  " << p << " " << u_x_t << endl;
    }

    return (u_x_t);
}
int main()
{
    system("chcp 1251");
    system("cls");
    
    cout << "Задача 4" << endl;
    cout << "Вариант 11" << endl;
   /*
    cout << "|---------------|----------------|----------------|----------------|----------------|\n";
    cout << "|     (h,tau)   |   (0.2,0.02)   |  (0.1,0.005)   | (0.05,0.00125) | (0.05,0.0052)  |\n";
    cout << "|---------------|----------------|----------------|----------------|----------------|\n";
    cout << "|---------------|----------------|----------------|----------------|----------------|\n";
   */
    double T = 0.1;
    int N;
    cout << "Введите число разбиений отрезка [0,1] -> ";
    cin >> N;
    int flag;
    cout << "Введите значение параметра \sigma для схемы с весами -> ";
    cin >> flag;
    cout << "Сетка значений точного решения:" << endl;
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    cout << "|      t/x      |       0       |       0.2     |      0.4      |      0.6      |       0.8     |       1       |\n";
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    for(int i = 0; i <=5; i++)
    {
        cout<<fixed<<setprecision(4)<< "|     "<<(T*i)/5<<"    |     " << accurate(N, 0, (T * i) / 5) << "    |     " << accurate(N, 0.2, (T * i) / 5) << "    |     " << accurate(N, 0.4, (T * i) / 5) << "    |     " << accurate(N, 0.6, (T * i) / 5) << "    |     " << accurate(N, 0.8, (T * i) / 5) << "    |     " << accurate(N, 1, (T * i) / 5) << "    |\n";
    }
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    cout << "Сетка значений решения, полученная по явной схеме с исп. ДРФ:" << endl;
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    cout << "|      t/x      |       0       |       0.2     |      0.4      |      0.6      |       0.8     |       1       |\n";
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    for (int i = 0; i <= 5; i++)
    {
        cout << fixed << setprecision(4) << "|     " << (T * i) / 5 << "    |     " << DFT(N, 0, (T * i) / 5) << "    |     " << DFT(N, 0.2, (T * i) / 5) << "    |     " << DFT(N, 0.4, (T * i) / 5) << "    |     " << DFT(N, 0.6, (T * i) / 5) << "    |     " << DFT(N, 0.8, (T * i) / 5) << "    |     " << DFT(N, 1, (T * i) / 5) << "    |\n";
    }
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    cout << "Сетка значений решения, полученная по схеме с весами с  исп. ДРФ:" << endl;    
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    cout << "|      t/x      |       0       |       0.2     |      0.4      |      0.6      |       0.8     |       1       |\n";
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    for (int i = 0; i <= 5; i++)
    {
        cout << fixed << setprecision(4) << "|     " << (T * i) / 5 << "    |     " << DFT_NET(flag, N, 0, (T * i) / 5) << "    |     " << DFT_NET(flag, N, 0.2, (T * i) / 5) << "    |     " << DFT_NET(flag, N, 0.4, (T * i) / 5) << "    |     " << DFT_NET(flag, N, 0.6, (T * i) / 5) << "    |     " << DFT_NET(flag, N, 0.8, (T * i) / 5) << "    |     " << DFT_NET(flag, N, 1, (T * i) / 5) << "    |\n";
    }
    cout << "|---------------|---------------|---------------|---------------|---------------|---------------|---------------|\n";
    cout << "Теперь подробнее про устойчивость:"<<endl;
    int flag2 = 1;
    int N2;
    cout << "Введите число разбиений отрезка [0,1] -> ";
    cin >> N2;
    cout << "Введите значение параметра \sigma для схемы с весами -> ";
    cin >> flag2;
    cout << DFT_NET_2(flag2, N2, 0.8, (T * 5) / 5);
    return 0;
}
