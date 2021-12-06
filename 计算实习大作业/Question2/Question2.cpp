// Question2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

constexpr auto N = 10;      // 矩阵维数
constexpr auto L = 100000;  // 迭代最大次数
constexpr auto epsilon = 1e-12; // 精度水平

double A[N][N] = { 0 };       // 定义实矩阵A
double Temp[N][N] = { 0 };    // 定义一个矩阵用于计算

void UpHessenberg(double T[N][N]);    // 拟上三角化
int Sgn(double d);     // 符号函数
void PrintAij(double T[N][N]);      // 在终端打印矩阵
void QRsolve(double A[N][N], double Mk[N][N], int m);       // QR方法
int DoubleStepQRsolve(double T[N][N], double Eigenvalue[N][2]);  // 双步位移QR方法
void SolveRoots(double a, double b, double c, double lamda[2][2]); // 求解一元二次方程的根
void GetMk(double T[N][N], double s, double t, double M_k[N][N], int m);   // 获取矩阵的M_k矩阵

// 符号函数
int Sgn(double d) {
    if (d < 0) return -1;
    else return 1;
}

// 在终端打印矩阵
void PrintAij(double T[N][N]) {
    cout << setiosflags(ios::scientific);   // 输出E型数
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            cout << setprecision(12) << T[i - 1][j - 1] << " ";
        }
        cout << endl;
    }
}

// 求解一元二次方程的根
void SolveRoots(double a, double b, double c, double lamda[2][2]) {
    // 一元二次方程的系数a,b,c
    // lamda[2][2]保存的方程的实根/复根
    double delta = b * b - 4 * a * c;
    if (delta >= 0) {// 方程有两个实根
        lamda[0][0] = (-b - sqrt(delta)) / (2 * a);
        lamda[0][1] = 0;
        lamda[1][0] = (-b + sqrt(delta)) / (2 * a);
        lamda[1][1] = 0;
    }
    else {
        lamda[0][0] = -b / (2 * a);
        lamda[0][1] = -sqrt(-delta) / (2 * a);
        lamda[1][0] = -b / (2 * a);
        lamda[1][1] = sqrt(-delta) / (2 * a);
    }
}

// 获取矩阵的M_k矩阵
void GetMk(double T[N][N], double s, double t, double M_k[N][N], int m) {
    int i, j;
    int k, h;
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            double sum = 0;
            for (k = 0; k < m; k++) {
                for (h = 0; h < m; h++) {
                    sum += T[i][k] * T[h][j];
                }
            }
            M_k[i][j] = sum - s * T[i][j];
            if (i == j) M_k[i][j] += t;
        }
    }
}

// QR方法
void QRsolve(double A[N][N], double Mk[N][N], int m) {
    int i, j, r;
    double dr, cr, hr, tr;
    double sum;
    double B[N][N], C[N][N];
    double ur[N], vr[N], pr[N], qr[N], wr[N];

    for (i = 1; i <= m; i++) {
        for (j = 1; j < +m; j++) {
            B[i - 1][j - 1] = Mk[i - 1][j - 1];
            C[i - 1][j - 1] = A[i - 1][j - 1];
        }
    }

    r = 1;
    while (r<=m-1) {
    Sentence_1:
        int flag = 0;
        for (i = r + 1; i <= m; i++) {
            if (B[i - 1][r - 1] == 0) flag = 0;
            else flag = 1;
        }
        if (flag == 0)
            goto Sentence_5;
        else goto Sentence_2;

    Sentence_2:
        dr = 0;
        for (i = r; i <= m; i++) dr += B[i - 1][r - 1] * B[i - 1][r - 1];
        dr = sqrt(dr);
        if (B[r - 1][r - 1] == 0) cr = dr;
        else cr = -Sgn(B[r - 1][r - 1]) * dr;
        hr = cr * cr - cr * B[r - 1][r - 1];
        for (i = 1; i <= m; i++) {
            if (i < r) ur[i - 1] = 0;
            else if (i == r) ur[i - 1] = B[r - 1][r - 1] - cr;
            else ur[i - 1] = B[i - 1][r - 1];
        }
        for (i = 1; i <= m; i++) {
            sum = 0;
            for (j = 1; j <= m; j++) sum += B[j - 1][i - 1] * ur[j - 1];
            vr[i - 1] = sum / hr;
        }
        for (i = 1; i <= m; i++)
            for (j = 1; j <= m; j++)
                B[i - 1][j - 1] = B[i - 1][j - 1] - ur[i - 1] * vr[j - 1];
        for (i = 1; i <= m; i++) {
            sum = 0;
            for (j = 1; j <= m; j++) sum += C[j - 1][i - 1] * ur[j - 1];
            pr[i - 1] = sum / hr;
        }
        for (i = 1; i <= m; i++) {
            sum = 0;
            for (j = 1; j <= m; j++) sum += C[i - 1][j - 1] * ur[j - 1];
            qr[i - 1] = sum / hr;
        }
        sum = 0;
        for (i = 1; i <= m; i++) sum += pr[i - 1] * ur[i - 1];
        tr = sum / hr;
        for (i = 1; i <= m; i++) wr[i - 1] = qr[i - 1] - tr * ur[i - 1];
        for (i = 1; i <= m; i++)
            for (j = 1; j <= m; j++)
                C[i - 1][j - 1] = C[i - 1][j - 1] - wr[i - 1] * ur[j - 1] - ur[i - 1] * pr[j - 1];
        goto Sentence_5;

    Sentence_5:
        r++;

    }

    for (i = 1; i <= m; i++)
        for (j = 1; j <= m; j++)
            A[i - 1][j - 1] = C[i - 1][j - 1];
}

// 拟上三角化
void UpHessenberg(double T[N][N]) {
    int i, j, r;
    double d_r, c_r, h_r, t_r, sum;
    double u_r[N], p_r[N], q_r[N], w_r[N];
    for (r = 1; r <= N - 2; r++) {
        int flTg = 0;   // 判断全为零的标志
        for (i = r + 2; i <= N; i++) {
            if (T[i - 1][r - 1] != 0) {
                flTg = 1;
                break;
            }
        }
        if (flTg == 0) continue;

        d_r = 0;    // d_r置零
        for (i = r + 1; i <= N; i++) d_r += T[i - 1][r - 1] * T[i - 1][r - 1];
        d_r = sqrt(d_r);
        if (T[r][r - 1] == 0) c_r = d_r;
        else c_r = -Sgn(T[r][r - 1]) * d_r;
        h_r = c_r * c_r - c_r * T[r][r - 1];
        for (i = 1; i <= N; i++) {
            if (i <= r) u_r[i - 1] = 0;
            else if (i == r + 1) u_r[i - 1] = T[r][r - 1] - c_r;
            else u_r[i - 1] = T[i - 1][r - 1];
        }
        for (i = 1; i <= N; i++) {
            sum = 0;
            for (j = 1; j <= N; j++) sum += T[j - 1][i - 1] * u_r[j - 1];
            p_r[i - 1] = sum / h_r;
        }
        for (i = 1; i <= N; i++) {
            sum = 0;
            for (j = 1; j <= N; j++) sum += T[i - 1][j - 1] * u_r[j - 1];
            q_r[i - 1] = sum / h_r;
        }
        sum = 0;
        for (i = 1; i <= N; i++) {
            sum += p_r[i - 1] * u_r[i - 1];
        }
        t_r = sum / h_r;
        for (i = 1; i <= N; i++) {
            w_r[i - 1] = q_r[i - 1] - t_r * u_r[i - 1];
        }
        // 求T(r+1)
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                T[i - 1][j - 1] = T[i - 1][j - 1] - w_r[i - 1] * u_r[j - 1] - u_r[i - 1] * p_r[j - 1];
            }
        }
    }
}

int DoubleStepQRsolve(double T[N][N], double Eigenvalue[N][2]) {

    int k = 0, m = N;
    int flag = 0;
    double M_k[N][N] = { 0 };
    double s, t, b, c;

    Sentence_3:
        if (fabs(T[m - 1][m - 2]) <= epsilon) {
            // 得到A的一个特征值
            Eigenvalue[flag][0] = T[m - 1][m - 1];
            flag++;
            m = m - 1;
            goto Sentence_4;
        }
        else goto Sentence_5;

    Sentence_4:
        if (m == 1) {
            Eigenvalue[flag][0] = T[0][0];
            flag++;
            goto Sentence_11;
        }
        else if (m == 0)
            goto Sentence_11;
        else if(m>1)
            goto Sentence_3;

    Sentence_5:
        b = -(T[m - 2][m - 2] + T[m - 1][m - 1]);
        c = T[m - 2][m - 2] * T[m - 1][m - 1] - T[m - 1][m - 2] * T[m - 2][m - 1];
        double lamda[2][2];
        SolveRoots(1, b, c, lamda);
        if (m == 2) {
            // 得到A的两个特征值
            Eigenvalue[flag][0] = lamda[0][0]; Eigenvalue[flag][1] = lamda[0][1];
            flag++;
            Eigenvalue[flag][0] = lamda[1][0]; Eigenvalue[flag][1] = lamda[1][1];
            flag++;
            goto Sentence_11;
        }
        else goto Sentence_7;

    Sentence_7:
        if (fabs(T[m - 2][m - 3]) <= epsilon) {
            // 得到A的两个特征值
            Eigenvalue[flag][0] = lamda[0][0]; Eigenvalue[flag][1] = lamda[0][1];
            flag++;
            Eigenvalue[flag][0] = lamda[1][0]; Eigenvalue[flag][1] = lamda[1][1];
            flag++;
            m = m - 2;
            goto Sentence_4;
        }
        else goto Sentence_8;

    Sentence_8:
        if (k == L) {
            cout << "计算终止，未得到A的全部特征值" << endl;
            return 0;
        }
        else goto Sentence_9;

    Sentence_9:
        s = T[m - 2][m - 2] + T[m - 1][m - 1];
        t = T[m - 2][m - 2] * T[m - 1][m - 1] - T[m - 1][m - 2] * T[m - 2][m - 1];
        // 获取A_k矩阵的M_k矩阵
        GetMk(T, s, t, M_k, m);
        // 对M_k作QR分解
        QRsolve(T, M_k, m);
        k = k + 1;
        goto Sentence_3;

    Sentence_11:
        cout << "A的全部特征值已计算完毕，停止计算" << endl;
        return 0;
    
}

int main()
{
    std::cout << "Hello World!\n";

    int i, j;

    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            if (i != j) A[i-1][j-1] = sin(0.5 * i + 0.2 * j);
            else A[i-1][j-1] = 1.52 * cos(i + 1.2 * j);
        }
    }
    cout << "实矩阵A为: " << endl;
    //PrintAij(A);

    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++)
            Temp[i - 1][j - 1] = A[i - 1][j - 1];
    }

    UpHessenberg(Temp);
    cout << "对A拟上三角化： " << endl;
    PrintAij(Temp);

    // 双步位移QR方法
    double Eigenvalue[N][2];
    for (i = 1; i <= N; i++)
        for (j = 1; j <= 2; j++)
            Eigenvalue[i - 1][j - 1] = 0;
    DoubleStepQRsolve(Temp, Eigenvalue);
    // 在终端打印特征值
    cout << setiosflags(ios::scientific);   // 输出E型数
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= 2; j++) {
            if (Eigenvalue[i - 1][0] != 0)
                cout << setprecision(12) << Eigenvalue[i - 1][j - 1] << " ";
        }
        cout << endl;
    }


    return 0;
}

