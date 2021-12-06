// Question1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 数值分析计算实习大作业，第一题

#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

#define N 501   // 矩阵阶数

double sigma = 1e-12;   // 精度水平=1e-12

double A[5][N] = {0}; // 定义对称矩阵A
double Tem[5][N] = { 0 };

void LUresolve() {// LU分解法
    int i, j, p, k, t;
    for (i = 3; i < 5; i++)
        Tem[i][0] = Tem[i][0] / Tem[2][0];
    for (k = 2; k <= N; k++) {
        if (k + 2 < N)
            p = k + 2;
        else p = N;
        for (j = k; j <= p; j++) {
            t = 1;
            if (k - 2 > t) t = k - 2;
            else t = t;
            if (j - 2 > t) t = j - 2;
            else t = t;
            double l = 0;
            for (; t < k; t++)
                l = l + Tem[k - t + 2][t - 1] * Tem[t - j + 2][j - 1];
            Tem[k - j + 2][j - 1] = Tem[k - j + 2][j - 1] - l;
        }
        if (k < N) {
            if (k + 2 < N) p = k + 2;
            else p = N;
            for (i = k + 1; i <= p; i++) {
                t = 1;
                if (i - 2 > t)t = i - 2;
                else t = t;
                if (k - 2 > t) t = k - 2;
                else t = t;
                double l = 0;
                for (; t < k; t++)
                    l = l + Tem[i - t + 2][t - 1] * Tem[t - k + 2][k - 1];
                Tem[i - k + 2][k - 1] = Tem[i - k + 2][k - 1] - l;
                Tem[i - k + 2][k - 1] = Tem[i - k + 2][k - 1] / Tem[2][k - 1];
            }
        }
    }
}

double PowerMethod() {// 幂法函数
    int i, j, p, s;
    double Temp, Mo;
    double u[N] = { 0 };
    double y[N] = { 0 };
    double beta1 = 0;
    double beta2 = 0;

    for (i = 0; i < N; i++)
        u[i] = { 1 };
    while (1) {// 开始迭代
        Temp = 0;
        for (i = 0; i < N; i++)
            Temp = Temp + u[i] * u[i];
        Mo = sqrt(Temp);    // 得到2范数
        // 归一化，得到向量y
        for (i = 0; i < N; i++)
            y[i] = u[i] / Mo;
        for (i = 0; i < N; i++)
            u[i] = 0;
        // 计算u_k
        for (i = 2; i < N + 2; i++) {
            if (i < 4)
                p = i;
            else p = 4;
            if (i < N - 1)
                s = i;
            else s = N - 1;
            for (j = i - p; j <= s; j++)
                u[i - 2] = u[i - 2] + Tem[i - j][j] * y[j];
        }
        // 计算本次迭代的Beta
        for (i = 0; i < N; i++)
            beta2 = beta2 + u[i] * y[i];
        if (fabs(beta2 - beta1) / fabs(beta2) <= sigma)
            return beta2;   // 达到精度要求，返回beta2
        else {// 转入下次迭代
            beta1 = beta2;
            beta2 = 0;
            //cout << "beta1=" << beta1 << endl;
        }
    }
}

double FPowerMethod() {// 反幂法函数
    int i, p, s;
    double Temp, Mo;
    double u[N] = { 0 };
    double y[N] = { 0 };
    double D[N] = { 0 };
    double beta1 = 0;
    double beta2 = 0;
    for (i = 0; i < N; i++)
        u[i] = 1;
    // 进行LU分解
    LUresolve();
    while (1) {// 开始迭代
        Temp = 0;
        for (i = 0; i < N; i++)
            Temp = Temp + u[i] * u[i];
        Mo = sqrt(Temp);    // 得到2范数
        // 归一化，得到向量y
        for (i = 0; i < N; i++)
            y[i] = u[i] / Mo;
        for (i = 0; i < N; i++)
            u[i] = 0;
        // LU分解求解u_k
        D[0] = y[0];
        for (i = 1; i < N; i++) {
            s = 1;
            if (i - 1 > s) s = i - 1;
            else s = s;
            double l = 0;
            for (; s <= i; s++)
                l = l + Tem[i - s + 3][s - 1] * D[s - 1];
            D[i] = y[i] - l;
        }
        u[N - 1] = D[N - 1] / Tem[2][N - 1];
        for (i = N - 2; i >= 0; i--) {
            if (i + 2 >= N) p = N - 1;
            else p = i + 2;
            double l = 0;
            for (s = i + 1; s <= p; s++)
                l = l + Tem[i - s + 2][s] * u[s];
            u[i] = (D[i] - l) / Tem[2][i];
        }
        for (i = 0; i < N; i++)
            beta2 = beta2 + y[i] * u[i];
        if(fabs(beta2-beta1)/fabs(beta2)<=sigma)
            return 1 / beta2;   // 达到精度要求，返回beta2
        else {// 转入下次迭代
            beta1 = beta2;
            beta2 = 0;
            //cout << "beta1=" << beta1 << endl;
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

    // 第一问
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            Tem[i][j] = A[i][j];
    }
    double Lam1 = PowerMethod();    // 幂法求最大特征值LamD
    //cout << "Lam1=" << Lam1 << endl;

    // 构造带原点位移的矩阵
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++) {
            if (i == 2)
                Tem[i][j] = A[i][j] - Lam1;
            else
                Tem[i][j] = A[i][j];
        }
    }
    double Lam2 = PowerMethod() + Lam1;
    //cout << "Lam2=" << Lam2 << endl;

    double LamMin, LamMax;
    if (Lam1 < Lam2) {
        LamMin = Lam1;
        LamMax = Lam2;
    }
    else {
        LamMin = Lam2;
        LamMax = Lam1;
    }

    cout << setiosflags(ios::scientific);
    cout << setprecision(12) << "Lam_1=" << LamMin << endl;
    cout <<setprecision(12) << "Lam_501=" << LamMax << endl;

    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            Tem[i][j] = A[i][j];
    }
    double Lam_s = FPowerMethod();
    cout <<setprecision(12) << "Lam_s=" << Lam_s << endl;

    //第二问
    for (k = 1; k < 40; k++) {
        double U_k = LamMin + k * (LamMax - LamMin) / 40;
        // 构造带位移的矩阵Tem
        for (i = 0; i < 5; i++) {
            for (j = 0; j < N; j++) {
                if (i == 2)
                    Tem[i][j] = A[i][j] - U_k;
                else
                    Tem[i][j] = A[i][j];
            }
        }
        double Lam_ik = FPowerMethod() + U_k;
        cout << "when k = " << k << "  " << setprecision(12) << "Lam_ik=" << Lam_ik << endl;
    }

    // 第三问
    cout << "cond(A)2 = " << setprecision(12) << fabs(Lam1 / Lam_s) << endl;

    // LU分解
    for (i = 0; i < 5; i++) {
        for (j = 0; j < N; j++)
            Tem[i][j] = A[i][j];
    }
    LUresolve();
    double det = 1;
    for (i = 0; i < N; i++) // 得到对角线元素
        det = det * Tem[2][i];
    cout << "detA = " << setprecision(12) << det << endl;

}

