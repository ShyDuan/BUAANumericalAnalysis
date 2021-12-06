// Quesiton2_new.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

constexpr auto N = 10;      // 矩阵维数
constexpr auto L = 100000;  // 迭代最大次数
constexpr auto epsilon = 1e-12; // 精度水平

void UpHessenberg(double T[N][N]);    // 拟上三角化
void CreatA(double T[N][N]);    // 构造矩阵A
int Sgn(double d);     // 符号函数
void PrintAij(double T[N][N]);      // 在终端打印矩阵
void PrintAi(double x[N]);      // 在终端打印向量
void QRsolve(double A[N][N], double Mk[N][N], int m);       // QR方法
int DoubleStepQRsolve(double T[N][N], double Eigenvalue[N][2]);  // 双步位移QR方法
void SolveRoots(double a, double b, double c, double lamda[2][2]); // 求解一元二次方程的根
void GetMk(double T[N][N], double s, double t, double M_k[N][N], int m);   // 获取矩阵的M_k矩阵

// 构造矩阵A
void CreatA(double A[N][N]) {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            if (i != j) A[i - 1][j - 1] = sin(0.5 * i + 0.2 * j);
            else A[i - 1][j - 1] = 1.52 * cos(i + 1.2 * j);
        }
    }
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

// 在终端打印向量
void PrintAi(double x[N]) {
    cout << "[";
    for (int i = 0; i < N; i++) {
        cout << x[i];
        if (i < N - 1) cout << ",";
    }
    cout << "]" << endl;
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
    int k;

    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
            M_k[i][j] = 0;

    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            double sum = 0;
            for (k = 0; k < m; k++)
                sum += T[i][k] * T[k][j];
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
        for (j = 1; j <= m; j++) {
            B[i - 1][j - 1] = Mk[i - 1][j - 1];
            C[i - 1][j - 1] = A[i - 1][j - 1];
        }
    }

    for (r = 1; r < m; r++) {
        int flag = 0;
        for (i = r + 1; i <= m; i++) {
            if (B[i - 1][r - 1] != 0) {
                flag = 1;
                break;
            }
            else flag = 0;
        }
        if (flag != 0) {
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
        }
    }

    for (i = 1; i <= m; i++)
        for (j = 1; j <= m; j++)
            A[i - 1][j - 1] = C[i - 1][j - 1];
}

// 双步位移QR方法
int DoubleStepQRsolve(double T[N][N], double Eigenvalue[N][2]) {
    int k = 1, m = N;
    double Mk[N][N];

    while (k < L) {
        k = k + 1;
        if (m == 0) return 0;
        else if (m == 1) {// 矩阵有一个特征值
            Eigenvalue[m - 1][0] = T[0][0];
            m = 0;
            continue;
        }
        else if (m == 2) {// 矩阵有两个特征值
            double s = -(T[m - 2][m - 2] + T[m - 1][m - 1]);
            double t = T[m - 2][m - 2] * T[m - 1][m - 1] - T[m - 1][m - 2] * T[m - 2][m - 1];
            double lamda[2][2];
            SolveRoots(1, s, t, lamda);
            // 得到A的两个特征值
            Eigenvalue[m - 2][0] = lamda[0][0]; Eigenvalue[m - 2][1] = lamda[0][1];
            Eigenvalue[m - 1][0] = lamda[1][0]; Eigenvalue[m - 1][1] = lamda[1][1];
            m = 0;
            continue;
        }
        else if (fabs(T[m - 1][m - 2]) <= epsilon) {
            Eigenvalue[m - 1][0] = T[m - 1][m - 1]; // 得到A的一个特征值
            m = m - 1;
            continue;
        }
        else {
            double s = -(T[m - 2][m - 2] + T[m - 1][m - 1]);
            double t = T[m - 2][m - 2] * T[m - 1][m - 1] - T[m - 1][m - 2] * T[m - 2][m - 1];

            if (fabs(T[m - 2][m - 3]) <= epsilon) {
                double lamda[2][2];
                SolveRoots(1, s, t, lamda);
                // 得到A的两个特征值
                Eigenvalue[m - 2][0] = lamda[0][0]; Eigenvalue[m - 2][1] = lamda[0][1];
                Eigenvalue[m - 1][0] = lamda[1][0]; Eigenvalue[m - 1][1] = lamda[1][1];
                m = m - 2;
                continue;
            }
            else {
                // 获取A_k矩阵的M_k矩阵
                GetMk(T, s, t, Mk, m);
                // 对M_k作QR分解
                QRsolve(T, Mk, m);
            }
        }
        // cout << "k=" << k << endl;
    }
    cout << "Error! 达到迭代最大次数，k=" << k << endl;
    return 0;
}

// Gauss消去法
void Gauss(double T[N][N], double ans[N]) {
    int k, i, j;
    double temp;

    // 消元过程
    for (k = 0; k < N - 1; k++) {
        int flag = k;
        temp = fabs(T[k][k]);
        for (j = k + 1; j < N; j++)
            if (fabs(T[j][k]) > temp) {
                flag = j;
                temp = fabs(T[j][k]);   // 找出主元
            }
        if (flag != k)
            for (j = k; j < N; j++)
                swap(T[k][j], T[flag][j]);  // 交换主元所在行全部元素
        for (i = k + 1; i < N; i++) {// 消元
            temp = T[i][k] / T[k][k];
            for (j = k + 1; j < N; j++)
                T[i][j] -= temp * T[k][j];
        }
    }
    //PrintAij(T);
    // 回代过程
    for (k = N - 1, ans[N - 1] = 1; k >= 1; k--) {
        for (j = k + 1, temp = 0; j <= N; j++)
            temp += T[k - 1][j - 1] * ans[j - 1];
        ans[k - 1] = -temp / T[k - 1][k - 1];
    }
    // 特征向量归一化
    temp = 0;
    for (i = 0; i < N; i++)
        temp += pow(ans[i], 2);
    temp = sqrt(temp);
    for (i = 0; i < N; i++)
        ans[i] = ans[i] / temp;
}


int main() {
    std::cout << "Hello World!\n";

    int i, j;

    double A[N][N] = { 0 };       // 定义实矩阵A
    CreatA(A);
    //cout << "实矩阵A为: " << endl;
    //PrintAij(A);

    UpHessenberg(A);
    cout << "对A拟上三角化： " << endl;
    PrintAij(A);

    // 双步位移QR方法
    double Eigenvalue[N][2];
    for (i = 1; i <= N; i++)
        for (j = 1; j <= 2; j++)
            Eigenvalue[i - 1][j - 1] = 0;
    DoubleStepQRsolve(A, Eigenvalue);

    // 在终端打印特征值
    cout << "矩阵A的特征值： " << endl;
    for (i = 1; i <= N; i++) {
        cout << setiosflags(ios::scientific);   // 输出E型数
        if (Eigenvalue[i - 1][1] != 0)
            cout << "λ = " << i << ", 特征值为 (" << setprecision(12) << Eigenvalue[i - 1][0] << ", "
            << Eigenvalue[i - 1][1] << ")" << endl;
        else
            cout << "λ = " << i << ", 特征值为 " << setprecision(12) << Eigenvalue[i - 1][0] << endl;
    }

    // 求解实特征值对应的特征向量
    double Temp[N][N];  // 系数矩阵
    double x[N];        // 特征向量
    for (int k = 1; k <= N; k++) {
        if (Eigenvalue[k - 1][1] == 0) {// 实特征值
            // 构造矩阵A
            CreatA(A);
            // 构造矩阵Temp=A-λI
            for (i = 0; i < N; i++)
                for (j = 0; j < N; j++) {
                    if (i == j) Temp[i][j] = A[i][j] - Eigenvalue[k - 1][0];
                    else Temp[i][j] = A[i][j];
                }
            // 高斯消去法求解特征向量
            cout << "λ = " << setprecision(12) << Eigenvalue[k - 1][0] << endl;
            cout << "对应的特征向量为：" << endl;
            Gauss(Temp, x);     // 列主元Guass消元法
            PrintAi(x);         // 打印特征向量
        }
    }
}
