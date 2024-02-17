#include <stdio.h>
#include <math.h>
#pragma warning(disable:4996)
#define m 3			
#define n 5			
double M = 1e6;
double A[m][n];		
double C[n];		
double b[m];		
double CB[m];		
double seta[m];		
double sigma[n];	
double x[n];		
int num[m];		
double Z = 0;	

void input();
void print();
int FindSwapInVar();		
int FindSwapOutVar(int a);	
void Iterate(int a, int b);	


int FindSwapInVar()
{
	int i, k = 0;
	int flag = 1;
	double min = 0;
	for (i = 0; i < n; i++)
	{
		if (sigma[i] < 0)
		{
			flag = 0; break;
		}
	}
	if (flag == 1)
		return -1;
	for (i = 0; i < n; i++)
	{
		if (sigma[i] < min)
		{
			min = sigma[i];
			k = i;
		}
	}
	return k;
}

int FindSwapOutVar(int a)
{
	int i, j;
	int flag = 1;
	int k = a;
	for (i = 0; i < m; i++)
	{
		if (A[i][k] > 0)
		{
			flag = 0; break;
		}
	}
	if (flag == 1)
	{
		printf("该线性规划问题有无界解、无最优解！\n");
		return -1;
	}
	for (i = 0; i < m; i++)
	{
		if (A[i][k] > 0)
			seta[i] = b[i] / A[i][k];
		else seta[i] = M;
	}
	double min = M;
	for (i = 0; i < m; i++)
	{
		if (min >= seta[i])
		{
			min = seta[i];
			j = i;
		}
	}
	num[j] = k + 1;
	CB[j] = C[k];
	return j;
}

void Iterate(int p, int q)
{
	int i, j, r, c;
	r = p;
	c = q;
	double temp1 = A[r][c];
	double temp2, temp3;

	b[r] /= temp1;
	for (j = 0; j < n; j++)
		A[r][j] /= temp1;
	for (i = 0; i < m; i++)
	{
		if (i != r)
			if (A[i][c] != 0)
			{
				temp2 = A[i][c];
				b[i] -= temp2 * b[r];
				for (j = 0; j < n; j++)
					A[i][j] -= temp2 * A[r][j];
			}
	}

	temp3 = sigma[c];
	for (i = 0; i < n; i++)
	{
		sigma[i] -= A[r][i] * temp3;
	}

}


void input()
{
	int i, j;
	printf("请按顺序输入方程组的系数矩阵A（共%d行%d列）:\n", m, n);
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			scanf("%lf", &A[i][j]);
	printf("请按顺序输入初始基变量所在列数的数字代码矩阵：\n");
	for (i = 0; i < m; i++)
		scanf("%d", &num[i]);
	printf("请按顺序输入约束条件的资源系数b矩阵：\n");
	for (i = 0; i < m; i++)
		scanf("%lf", &b[i]);
	printf("请按顺序输入目标函数的价值系数C矩阵：\n");
	for (i = 0; i < n; i++)
		scanf("%lf", &C[i]);
	for (i = 0; i < n; i++)
		sigma[i] = C[i];
	for (i = 0; i < m; i++)
		CB[i] = C[num[i] - 1];
}

void print()
{
	int i, j;
	printf("\n--------------------------------------------------------------------------\n");
	for (i = 0; i < m; i++)
	{
		printf("%8.2f\tX(%d) %8.2f ", CB[i], num[i], b[i]);
		for (j = 0; j < n; j++)
			printf("%8.2f ", A[i][j]);
		if (i != m - 1)
			printf("\n");
	}
	printf("\n--------------------------------------------------------------------------\n");
	printf("\t\t σ          ");
	for (i = 0; i < n; i++)
		printf(" %8.2f", sigma[i]);
	printf("\n--------------------------------------------------------------------------\n");
}
int main()
{
	int i, j, k = 1;
	int p, q;
	input();
	printf("\n--------------------------------------------------------------------------\n");
	printf("\tCB\tXB\tb\t");
	for (i = 0; i < n; i++)
		printf(" X(%d)\t", i + 1);
	for (i = 0; i < n; i++)
		x[i] = 0;
	while (1)
	{
		q = FindSwapInVar();
		if (q == -1)
		{
			print();
			printf("\n恭喜您得到最优解！\n");
			printf("最有解为：");
			for (j = 0; j < m; j++)
				x[num[j] - 1] = b[j];
			for (i = 0; i < n; i++)
			{
				printf("x%d=%.2f ", i + 1, x[i]);
				Z = Z + x[i] * C[i];
			}
			printf("\n最优值为：Z* = %.2f", Z);
			break;
		}
		print();
		p = FindSwapOutVar(q);
		printf("\n进行第%d次迭代,迭代元位置为(%d,%d)\n", k++, p + 1, q + 1);
		if (q == -1) break;
		Iterate(p, q);
	}
}


