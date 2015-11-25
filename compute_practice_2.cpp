// compute_practice_2.cpp 
//带双步位移的QR方法

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double e = 1e-12;      //迭代的精度水平
const int Max = 1000;        //迭代的最大次数

//定义复数结构体
typedef struct{
    double Re;
    double Im;
}ComplexNumber;

//初始化矩阵
void initMat(double **& a)
{
    int i, j;
    a = (double **)malloc(10 * sizeof(double *));
    for (i = 0; i<10; i++)
        a[i] = (double *)malloc(11 * sizeof(double));
    for (i = 0; i<10; i++) {
        for (j = 0; j<10; j++) {
            if (i != j)
                a[i][j] = sin(0.5*(i + 1) + 0.2*(j + 1));
            else
                a[i][j] = 1.52*cos(i + 1 + 1.2*(j + 1));
        }
    }
}

//将矩阵中小于e的值全部赋值为0
void zeroMat(double **a)
{
    int i, j;
    for (i = 0; i<10; i++) {
        for (j = 0; j<10; j++) {
            if (a[i][j]<e&&a[i][j]>-e)
                a[i][j] = 0;
        }
    }
}

//输出矩阵
int printMat(double **a)
{
    int i, j;
    for (i = 0; i<10; i++)
    {
        for (j = 0; j<10; j++)
            printf("%.12e ", a[i][j]);
        printf("\n");
    }
    printf("\n\n\n");
    return 0;
}

//将矩阵拟上三角化
void hessenbergMat(double **a)       
{
    int r, i, j;
    double c, d, h, t, u[10], p[10], q[10], w[10];
    for (r = 0; r<8; r++) {
        zeroMat(a);
        c = 0; d = 0; h = 0;
        for (i = r + 2; i<10; i++)
            d += a[i][r] * a[i][r];
        if (d == 0)
            continue;
        else {
            d += a[r + 1][r] * a[r + 1][r];
            d = pow(d, 0.5);
            if (a[r + 1][r] != 0)
                c = -fabs(a[r + 1][r]) / a[r + 1][r] * d;
            else 
                c = d;
            h = c*c - c*a[r + 1][r];
            for (i = 0; i<r + 1; i++)
                u[i] = 0;
            u[r + 1] = a[r + 1][r] - c;
            for (i = r + 2; i<10; i++) 
                u[i] = a[i][r];
            for (i = 0; i<10; i++) {
                p[i] = 0; q[i] = 0;
                for (j = 0; j<10; j++) {
                    p[i] += a[j][i] * u[j] / h;
                    q[i] += a[i][j] * u[j] / h;
                }
            }
            t = 0;
            for (i = 0; i<10; i++)
                t += p[i] * u[i] / h;
            for (i = 0; i<10; i++)
                w[i] = q[i] - t*u[i];
            for (i = 0; i<10; i++) {
                for (j = 0; j<10; j++)
                    a[i][j] -= (w[i] * u[j] + u[i] * p[j]);
            }
        }
    }
    printf("A(n-1)\n");
    zeroMat(a);
    printMat(a);
}

//矩阵乘法
void muiltiplyMat(double **a, double **b, double **c, int m)
{    
    for (int i = 0; i<m; i++) {
        for (int j = 0; j<m; j++) {
            c[i][j] = 0;
            for (int k = 0; k<m; k++)
                c[i][j] += a[i][k] * b[k][j];
        }
    }
}

//QR分解中的迭代运算
void iterate(double **M, double **a, int m)
{
    int r, i, j;
    double c, d, h, t, u[10], v[10], p[10], q[10], w[10];
    for (r = 0; r<m - 1; r++) {
        zeroMat(M);
        c = 0; d = 0; h = 0;
        for (i = r + 1; i<m; i++)
            d += M[i][r] * M[i][r];
        if (fabs(d) == 0)
            continue;
        else {
            d += M[r][r] * M[r][r];
            d = pow(d, 0.5);
            if (M[r][r] != 0)
                c = -fabs(M[r][r]) / M[r][r] * d;
            else 
                c = d;
            h = c*c - c*M[r][r];
            for (i = 0; i<r; i++)u[i] = 0;
            u[r] = M[r][r] - c;
            for (i = r + 1; i<m; i++) 
                u[i] = M[i][r];
            for (i = 0; i<m; i++) {
                v[i] = 0;
                for (j = 0; j<m; j++)
                    v[i] += M[j][i] * u[j] / h;
            }
            for (i = 0; i<m; i++) {
                p[i] = 0; q[i] = 0;
                for (j = 0; j<m; j++) {
                    M[i][j] -= u[i] * v[j];
                    p[i] += a[j][i] * u[j] / h;
                    q[i] += a[i][j] * u[j] / h;
                }
            }
            t = 0;
            for (i = 0; i<m; i++)
                t += p[i] * u[i] / h;
            for (i = 0; i<m; i++)
                w[i] = q[i] - t*u[i];
            for (i = 0; i<m; i++) {
                for (j = 0; j<m; j++) {
                    a[i][j] -= (w[i] * u[j] + u[i] * p[j]);
                }
            }
        }
    }
}

//Gauss消元法中的选主元
void maxline(double **ptr, int n, int k)
{
    double c;
    int i, M;
    M = k;
    for (i = k; i<n; i++) {
        if (fabs(ptr[i][k])>fabs(ptr[M][k]))
            M = i;
    }
    if (M>k) {
        for (i = k; i<n + 1; i++) {
            c = ptr[k][i];
            ptr[k][i] = ptr[M][i];
            ptr[M][i] = c;
        }
    }
}

//Gauss消元法
void gauss(double lambda)
{
    double **a;
    double *X;
    double m, sigma, t;
    int i, j, k, n = 10;
    X = (double *)malloc(10 * sizeof(double));
    initMat(a);

    for (i = 0; i<n; i++) {
        a[i][i] -= lambda;
        a[i][10] = 0;
    }

    for (k = 0; k<n - 1; k++) {
        maxline(a, n, k);
        for (i = k + 1; i<n; i++) {
            m = a[i][k] / a[k][k];
            for (j = k; j<n + 1; j++)
                a[i][j] = a[i][j] - m*a[k][j];
        }
    }

    X[n - 1] = 1;    
    t = 0;
    for (k = n - 2; k >= 0; k--) {
        sigma = 0;
        for (j = k + 1; j<n; j++)
            sigma += a[k][j] * X[j];
        X[k] = (a[k][n] - sigma) / a[k][k];
    }
    for (i = 0; i<n; i++)
        t += X[i] * X[i];
    t = sqrt(t);
    printf("eigenvector = ( ");
    for (i = 0; i<n; i++)
        printf("%.12e ", X[i]/t);
    printf(")\n");
    for (i = 0; i<10; i++)
        free(a[i]);
    free(a);
    free(X);
}

double inline sgn(double n)
{
    if (n > 0) {
        return 1;
    }
    else if (n < 0) {
        return -1;
    }
    else
        return 0;
}

//求Q、R和RQ
void QR_and_RQ(double **a)
{
    int i, j, r;
    double c, d, h, w[10], p[10], u[10];
    double **Q, **R, **RQ;
    Q = (double **)malloc(10 * sizeof(double *));
    for (i = 0; i < 10; i++)
        Q[i] = (double *)malloc(11 * sizeof(double));
    R = (double **)malloc(10 * sizeof(double *));
    for (i = 0; i < 10; i++)
        R[i] = (double *)malloc(11 * sizeof(double));
    RQ = (double **)malloc(10 * sizeof(double *));
    for (i = 0; i<10; i++)
        RQ[i] = (double *)malloc(11 * sizeof(double));
    for(i = 0; i < 10; i++) {
        for(j = 0; j < 10; j++){
            if(i == j)
                Q[i][i] = 1;
            else
                Q[i][j] = 0;
        }
    }
    for(i = 0; i < 10; i++) {
        for(j = 0; j < 10; j++)
            R[i][j] = a[i][j];
    }
    zeroMat(R);
    for(r = 0; r < 9; r++){
        d = 0;
        for(i = r + 1; i < 10; i++)
            d += R[i][r] * R[i][r];
        if(fabs(d) == 0)
            continue;
        else {
            d += R[r][r] * R[r][r];
            d = sqrt(d);
            if (R[r][r] == 0)
                c = d;
            else
                c = -sgn(R[r][r]) * d;
            h = c * c - c * R[r][r];
            for(i = 0; i < r; i++)
                u[i] = 0;
            u[r] = R[r][r] - c;
            for(i = r + 1; i < 10; i++)
                u[i] = R[i][r];
            for(i = 0; i < 10; i++) {
                w[i] = 0;
                for(j = 0; j < 10; j++)
                    w[i] += Q[i][j] * u[j];
            }
            for(i = 0; i < 10; i++) {
                for(j = 0; j < 10; j++)
                    Q[i][j] -= w[i] * u[j] / h;
            }
            for(i = 0; i < 10; i++){
                p[i] = 0;
                for(j = 0; j < 10; j ++ )
                    p[i] += R[j][i] * u[j] / h;
            }
            for(i = 0; i < 10; i++){
                for(j = 0; j < 10; j++)
                    R[i][j] -= u[i] * p[j];
            }
        }

    }
    zeroMat(R);
    zeroMat(Q);
    zeroMat(RQ);
    printf("Q:\n");
    printMat(Q);
    printf("R:\n");
    printMat(R);
    muiltiplyMat(R, Q, RQ, 10);
    printf("RQ:\n");
    printMat(RQ);
    for (i = 0; i < 10; i++)
        free(Q[i]);
    free(Q);
    for (i = 0; i < 10; i++)
        free(R[i]);
    free(R);
    for (i = 0; i<10; i++)
        free(RQ[i]);
    free(RQ);
} 

//带双布位移的QR分解方法
void QRmethod(double **a)
{
    int k, m, i, j, r;
    double s, t, det;
    double **M;
    ComplexNumber L[10];    
    M = (double**)malloc(10 * sizeof(double *));
    for (k = 0; k<10; k++)
        M[k] = (double*)malloc(10 * sizeof(double));

    m = 9; r = 0;
    for (k = 0; k<Max; k++)
    {
        if (m == 0) {
            L[r].Re = a[m][m];
            L[r].Im = 0; break;
        }
        else if (m < 0) {
            break;
        }
        if (fabs(a[m][m - 1])<e) {
            L[r].Re = a[m][m];
            L[r].Im = 0;
            m--; r++;
        }
        else {
            if (m == 1) {
                det = (a[m][m] + a[m - 1][m - 1])*(a[m][m] + a[m - 1][m - 1]) - 4 * (a[m][m] * a[m - 1][m - 1] - a[m - 1][m] * a[m][m - 1]);
                if (det > 0) {
                    L[r].Re = (a[m][m] + a[m - 1][m - 1]) / 2 + sqrt(det) / 2;
                    L[r].Im = 0;
                    L[r + 1].Re = (a[m][m] + a[m - 1][m - 1]) / 2 - sqrt(det) / 2;
                    L[r + 1].Im = 0;
                }
                else {
                    L[r].Re = (a[m][m] + a[m - 1][m - 1]) / 2;
                    L[r].Im = sqrt(-det) / 2;
                    L[r + 1].Re = (a[m][m] + a[m - 1][m - 1]) / 2;
                    L[r + 1].Im = -sqrt(-det) / 2;
                }
                m -= 2;
                r += 2;
                continue;
            }
            else if (fabs(a[m - 1][m - 2])<e) {
                det = (a[m][m] + a[m - 1][m - 1])*(a[m][m] + a[m - 1][m - 1]) - 4 * (a[m][m] * a[m - 1][m - 1] - a[m - 1][m] * a[m][m - 1]);
                if (det>0) {
                    L[r].Re = (a[m][m] + a[m - 1][m - 1]) / 2 + sqrt(det) / 2;
                    L[r].Im = 0;
                    L[r + 1].Re = (a[m][m] + a[m - 1][m - 1]) / 2 - sqrt(det) / 2;
                    L[r + 1].Im = 0;
                }
                else {
                    L[r].Re = (a[m][m] + a[m - 1][m - 1]) / 2;
                    L[r].Im = sqrt(-det) / 2;
                    L[r + 1].Re = (a[m][m] + a[m - 1][m - 1]) / 2;
                    L[r + 1].Im = -sqrt(-det) / 2;
                }
                m -= 2;
                r += 2; 
                continue;
            }
            else {
                s = a[m - 1][m - 1] + a[m][m];
                t = a[m - 1][m - 1] * a[m][m] - a[m][m - 1] * a[m - 1][m];
                muiltiplyMat(a, a, M, m + 1);
                for (i = 0; i<10; i++)
                {
                    for (j = 0; j<10; j++)
                        M[i][j] -= s*a[i][j];
                    M[i][i] += t;
                }
                iterate(M, a, m + 1);
                zeroMat(a);
            }
        }
    }
    zeroMat(a);
    printf("after QR method\n");
    printMat(a);
    
    for (int r = 0; r<10; r++) {
        printf("\n");
        if (L[r].Im == 0) {
            printf("lambda[%d] = (%.12e + i*%.12e)\n", r + 1, L[r].Re, L[r].Im);
            gauss(L[r].Re);
        }
        else {
            printf("lambda[%d] = (%.12e + i*%.12e)\n", r + 1, L[r].Re, L[r].Im);
        }
    }
	for (int i = 0; i < 10; i++) {
		free(M[i]);
	}
	free(M);
}

int main()
{
    double **a;
    initMat(a);    
    hessenbergMat(a);    
    QR_and_RQ(a);
    QRmethod(a);
    return 0;
}
