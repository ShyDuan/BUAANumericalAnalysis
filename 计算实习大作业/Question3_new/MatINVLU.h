#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

//LUP�ֽ�
void LUP_Descomposition(vector<vector<double> >&A, vector<vector<double> >& L, vector<vector<double> >& U, int P[N])
{
    int row = 0;
    for (int i = 0; i < N; i++)
    {
        P[i] = i;
    }
    for (int i = 0; i < N - 1; i++)
    {
        double p = 0.000000;
        for (int j = i; j < N; j++)
        {
            if (abs(A[j * N + i]) > p)
            {
                p = abs(A[j * N + i]);
                row = j;
            }
        }
        if (0 == p)
        {
            cout << "�������죬�޷�������" << endl;
            return;
        }

        //����P[i]��P[row]
        int tmp = P[i];
        P[i] = P[row];
        P[row] = tmp;

        double tmp2 = 0.000000;
        for (int j = 0; j < N; j++)
        {
            //����A[i][j]�� A[row][j]
            tmp2 = A[i * N + j];
            A[i * N + j] = A[row * N + j];
            A[row * N + j] = tmp2;
        }

        //����ͬLU�ֽ�
        double u = A[i * N + i], l = 0.0;
        for (int j = i + 1; j < N; j++)
        {
            l = A[j * N + i] / u;
            A[j * N + i] = l;
            for (int k = i + 1; k < N; k++)
            {
                A[j * N + k] = A[j * N + k] - A[i * N + k] * l;
            }
        }

    }

    //����L��U
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                L[i * N + j] = A[i * N + j];
            }
            else
            {
                L[i * N + j] = 1;
            }
        }
        for (int k = i; k < N; k++)
        {
            U[i * N + k] = A[i * N + k];
        }
    }

}

//LUP��ⷽ��
vector<double> LUP_Solve(vector<vector<double> >L, vector<vector<double> >U, vector<int>P, vector<double>b)
{
    int N = L.size();
    vector<double>x(N);
    vector<double>y(N);

    //�����滻
    for (int i = 0; i < N; i++)
    {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++)
        {
            y[i] = y[i] - L[i * N][j] * y[j];
        }
    }
    //�����滻
    for (int i = N - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = N - 1; j > i; j--)
        {
            x[i] = x[i] - U[i * N][j] * x[j];
        }
        x[i] /= U[i * N][i];
    }
    return x;
}

/*****************����ԭ��ת��BEGIN********************/

/* ��� */
int getNext(int i, int m, int n)
{
    return (i % n) * m + i / n;
}

/* ǰ�� */
int getPre(int i, int m, int n)
{
    return (i % m) * n + i / m;
}

/* �������±�iΪ���Ļ� */
void movedata(double* mtx, int i, int m, int n)
{
    double temp = mtx[i]; // �ݴ�
    int cur = i;    // ��ǰ�±�
    int pre = getPre(cur, m, n);
    while (pre != i)
    {
        mtx[cur] = mtx[pre];
        cur = pre;
        pre = getPre(cur, m, n);
    }
    mtx[cur] = temp;
}

/* ת�ã���ѭ���������л� */
void transpose(double* mtx, int m, int n)
{
    for (int i = 0; i < m * n; ++i)
    {
        int next = getNext(i, m, n);
        while (next > i) // �����ں��С��i˵���ظ�,�Ͳ�������ȥ�ˣ�ֻ�в��ظ�ʱ����whileѭ����
            next = getNext(next, m, n);
        if (next == i)  // ����ǰ��
            movedata(mtx, i, m, n);
    }
}
/*****************����ԭ��ת��END********************/

//LUP����(��ÿ��b����ĸ���x������װ)
double* InvLU(vector<vector<double> > A, vector<vector<double> >& inv_A)
{
    vector<vector<double> >A_mirror(A); //��������A�ĸ�����ע�ⲻ��ֱ����A���㣬��ΪLUP�ֽ��㷨�ѽ���ı�
    vector<double>inv_A_each;   //������ĸ���
    vector<double>b;
    //double* b = new double[N]();//b��ΪB����о������
    int N = A.size();

    for (int i = 0; i < A.size(); i++)
    {
        vector<vector<double> >L(N, vector<double>(N));
        vector<vector<double> >U(N, vector<double>(N));
        vector<int>P(N);
        //double* L = new double[N * N]();
        //double* U = new double[N * N]();
        //int* P = new int[N]();

        //���쵥λ���ÿһ��
        for (int i = 0; i < N; i++)
        {
            b[i] = 0;
        }
        b[i] = 1;

        //ÿ�ζ���Ҫ���½�A����һ��
        for (int i = 0; i < N * N; i++)
        {
            A_mirror[i] = A[i];
        }

        LUP_Descomposition(A_mirror, L, U, P);

        inv_A_each = LUP_Solve(L, U, P, b);
        memcpy(inv_A + i * N, inv_A_each, N * sizeof(double));//������ƴ������
    }
    transpose(inv_A, N, N);//�������ڸ���ÿ��b�����x���д洢�������ת��

    return inv_A;
}
