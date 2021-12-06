// Question1_new.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 编程环境：VS2019

#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

constexpr auto N = 501;   // 矩阵阶数;;

double epsilon = 1e-12;   // 精度水平=1e-12

double A[5][N] = { 0 }; // 定义对称矩阵A
double Temp[5][N] = { 0 };

// 函数声明
double PowerMethod(double T[][N]);
double AntiPowerMethod(double T[][N]);
void LUresolve(double T[][N]);
int GetMin(int i, int j);
int GetMax(int i, int j);

int GetMin(int i, int j) {
    if (i <= j) return i;
    else return j;
}

int GetMax(int i, int j) {
    if (i > j) return i;
    else return j;
}

// 幂法的实现
double PowerMethod(double T[][N]) {
    int i, j, m, n;
    double u[N] = { 0 };    // 迭代向量
    double y[N] = { 0 };    // 归一化迭代向量
    double beta_last = 0;   // k-1时的β
    double beta_now = 0;    // k时的β
    for (i = 0; i < N; i++)
        u[i] = 1;   // 初始化迭代向量
    while (1) {// 开始迭代
        double Sum = 0;
        for (i = 0; i < N; i++)
            Sum += u[i] * u[i];
        double Mo = sqrt(Sum);  // 得到u(k-1)的模
        // 归一化迭代向量，得到y(k-1)
        for (i = 0; i < N; i++)
            y[i] = u[i] / Mo;
        for (i = 0; i < N; i++)
            u[i] = 0;

        // 计算u(k)：A与y(k-1)的乘积
        for (i = 2; i < N + 2; i++) {
            m = GetMin(i, 4);
            n = GetMin(i, N);
            for (j = i - m; j < n; j++)
                u[i - 2] = u[i - 2] + T[i - j][j] * y[j];
        }
        // 计算k-1时的β：yT(k-1)和u(k)的乘积
        for (i = 0; i < N; i++)
            beta_now = beta_now + y[i] * u[i];
        if (fabs(beta_now - beta_last) / fabs(beta_now) <= epsilon)
            return beta_now;    // 迭代结束，得到主特征值
        else {
            beta_last = beta_now;
            beta_now = 0;
            // cout << "beta_last = " << beta_last << endl;
        }
    }
}

// 反幂法的实现
double AntiPowerMethod(double T[][N]) {
    int i, j, k, m, n;
    double u[N] = { 0 };    // 迭代向量
    double y[N] = { 0 };    // 归一化的迭代向量
    double Temp_1[N] = { 0 };
    double beta_last = 0;
    double beta_now = 0;

    for (i = 0; i < N; i++)
        u[i] = 1;   // 初始化迭代向量
    // 对T作LU分解
    LUresolve(T);
    while (1) {
        double Sum = 0;
        for (i = 0; i < N; i++)
            Sum += u[i] * u[i];
        double Mo = sqrt(Sum);  // 得到u(k-1)的模
        // cout << "Mo = " << Mo << endl;
        // 归一化迭代向量，得到y(k-1)
        for (i = 0; i < N; i++)
            y[i] = u[i] / Mo;
        for (i = 0; i < N; i++)
            u[i] = 0;
        // 求解Ly=b
        Temp_1[0] = y[0];
        for (i = 1; i < N; i++) {
            double sum = 0;
            for (m = GetMax(1, i - 1); m <= i; m++)
                sum += T[i - m + 3][m - 1] * Temp_1[m - 1];
            Temp_1[i] = y[i] - sum;
        }
        // 求解Ux=y
        u[N - 1] = Temp_1[N - 1] / T[2][N - 1];
        for (i = N - 2; i >= 0; i--) {
            if (i + 2 >= N) n = N - 1;
            else n = i + 2;
            double sum = 0;
            for (j = i + 1; j <= n; j++)
                sum += T[i - j + 2][j] * u[j];
            u[i] = (Temp_1[i] - sum) / T[2][i];
        }
        // 计算k-1时的β：yT(k-1)和u(k)的乘积
        for (i = 0; i < N; i++)
            beta_now = beta_now + y[i] * u[i];
        if (fabs(beta_now - beta_last) / fabs(beta_now) <= epsilon)
            return 1 / beta_now;    // 迭代结束，得到主特征值
        else {
            beta_last = beta_now;
            beta_now = 0;
            // cout << "beta_last = " << beta_last << endl;
        }

    }
}

// LU分解的实现
void LUresolve(double T[][N]) {
    int i, j, p, k, t;
    for (i = 3; i < 5; i++)
        T[i][0] = T[i][0] / T[2][0];
    for (k = 2; k <= N; k++) {
        p = GetMin(k + 2, N);
        for (j = k; j <= p; j++) {// 得到U
            double sum = 0;
            for (t = GetMax(1, GetMax(k - 2, j - 2)); t < k; t++)
                sum += T[k - t + 2][t - 1] * T[t - j + 2][j - 1];
            T[k - j + 2][j - 1] = T[k - j + 2][j - 1] - sum;
        }
        if (k < N) {
            p = GetMin(k + 2, N);
            for (i = k + 1; i <= p; i++) {
                double sum = 0;
                for (t = GetMax(1, GetMax(i - 2, k - 2)); t < k; t++)
                    sum += T[i - t + 2][t - 1] * T[t - k + 2][k - 1];
                T[i - k + 2][k - 1] = T[i - k + 2][k - 1] - sum;
                T[i - k + 2][k - 1] = T[i - k + 2][k - 1] / T[2][k - 1];
            }
        }
    }
}

int main()
{
    cout << "Hello World!\n";
    int i, j, k;
    for (i = 0; i < N; i++) {
        // 构造对称矩阵A
        k = i + 1;
        A[0][i] = -0.064;
        A[1][i] = 0.16;
        A[2][i] = (1.64 - 0.024 * k) * sin(0.2 * k) - 0.64 * exp(0.1 / k);
        A[3][i] = 0.16;
        A[4][i] = -0.064;
    }
    A[0][0] = A[0][1] = A[1][0] = 0;
    A[3][N - 1] = A[4][N - 1] = 0;
    A[4][N - 2] = 0;

    // 复制矩阵A，用于求解
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            Temp[i][j] = A[i][j];
    }
    
    // 幂法求解主特征值
    double Lamd1 = PowerMethod(Temp);
    
    // 构造带原点平移的矩阵Temp
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++) {
            if (i == 2) Temp[i][j] = A[i][j] - Lamd1;
            else Temp[i][j] = A[i][j];
        }
    }
    // Lamd1模值最大，为最大值最小值之一
    // Lamd2离Lamd1最远，故也为最大值最小值之一
    double Lamd2 = PowerMethod(Temp) + Lamd1;

    cout << setiosflags(ios::scientific);   // 输出E型数
    // cout << setprecision(12) << "Lamd1 = " << Lamd1 << endl;
    // cout << setprecision(12) << "Lamd2 = " << Lamd2 << endl;
    
    double Lamd_1, Lamd_501;
    if (Lamd1 >= Lamd2) {// 比较特征值大小
        Lamd_1 = Lamd2;
        Lamd_501 = Lamd1;
    }
    else {
        Lamd_1 = Lamd1;
        Lamd_501 = Lamd2;
    }
    cout << setprecision(12) << "Lamd_1 = " << Lamd_1 << endl;
    cout << setprecision(12) << "Lamd_501 = " << Lamd_501 << endl;

    // 复制矩阵A，用于求解
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            Temp[i][j] = A[i][j];
    }

    // 反幂法求解模最小的特征值
    double Lamd_s = AntiPowerMethod(Temp);
    cout << setprecision(12) << "Lamd_s = " << Lamd_s << endl;

    // 第二问
    for (k = 1; k < 40; k++) {
        double U_k = Lamd_1 + k * (Lamd_501 - Lamd_1) / 40;
        // 构造带位移的矩阵Tem
        for (i = 0; i < 5; i++) {
            for (j = 0; j < N; j++) {
                if (i == 2)
                    Temp[i][j] = A[i][j] - U_k;
                else
                    Temp[i][j] = A[i][j];
            }
        }
        double Lam_ik = AntiPowerMethod(Temp) + U_k;
        cout << "when k = " << k << "  " << setprecision(12) << "Lam_ik=" << Lam_ik << endl;
    }

    // 第三问
    cout << "cond(A)2 = " << setprecision(12) << fabs(Lamd1 / Lamd_s) << endl;

    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            Temp[i][j] = A[i][j];
    }
    LUresolve(Temp);
    double det = 1;
    for (i = 0; i < N; i++) det = det * Temp[2][i];
    cout << "detA = " << setprecision(12) << det << endl;
    return 0;

}
