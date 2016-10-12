#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES // for C
#include <math.h>
#define eps 0.00001;

const int N = 6;
double *X;
double *resultat;

void printmas(double *x, int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		printf("%f\n", x[i]);
	}
}

double kepler(const double *x, double *result)
//where x is array of [p, e, omega, i, Om, tau]
//result is array x, y, z, vx,vy, vz
{
	const int N = 100000;
	double u;//pericentr
	double r;//range
	double v;
	double temp[6];
	double a;
	double nu = 398600;//km^3/sec
	a = x[0] / (1 - x[1]);
	if (x[1] < 1)
	{
		int i;
		double f_n1;
		double f_n2;
		double E;
		double M;//&
		double t;//?
		M = sqrt(nu / (a * a * a))*(t - x[5]);
		f_n1 = M;
		for(i = 0; i < N; i++)
		{
			
			f_n2 = M + x[1]*sin(f_n1);
			f_n1 = f_n2;
			if (fabs(f_n2 - f_n1) < eps) break;
		}
		E = f_n2;
		v = 2 * arctg(sqrt((1 + x[1]) / (1 - x[1])) * sin(E) / (1 + cos(E)));
	}
	else if (x[1] == 1)
	{ 
		double M;//&
		M = sqrt(nu / (a * a * a)) * (t - x[5]);
		double x_ = 0.5 * pow((12 * M + 4 * sqrt((9 * M * M) + 4)), 1.0/3.0);
		v = 2 * arctg(x_ + (1.0/x_));
	}
	else
	{
		int i;
		double f_n1;
		double f_n2;
		double H;
		double M;//?
		double t;//?
		M = sqrt(nu / pow(x[0]/(x[1] * x[1] - 1),3))*(t - x[5]);
		f_n1 = M;
		for(i = 0; i < N; i++)
		{
			
			f_n2 = -M + x[1] * sinh(f_n1);
			f_n1 = f_n2;
			if (fabs(f_n2 - f_n1) < eps) break;
		}
		H = f_n2;
		v = 2 * arctg(sqrt((1 + x[1]) / (-1 + x[1])) * tanh(H / 2);
	}
	u = x[2] + v;
	r = x[0] / (1 + cos(v));
	temp[0] = (cos(u) * cos(x[4]) - sin(u) * sin(x[4]) * cos(x[3]));
	temp[1] = (cos(u) * sin(x[4]) + sin(u) * cos(x[4]) * cos(x[3]));
	temp[2] = sin(u) * sin(x[3]);
	result[0] = r * temp[0];//x
	result[1] = r * temp[1];//y
	result[2] = r * temp[2];//z
	temp[3] = (-sin(u) * cos(x[4]) - cos(u) * sin(x[4]) * cos(x[3]));
	temp[4] = (-sin(u) * sin(x[4]) + cos(u) * cos(x[4]) * cos(x[3]));
	temp[5] = (x[1] * sin(v) * sin(u) * sin(x[4]) + (1 + x[1] * cos(v)) * cos(u) * cos(x[4]));
	result[3] = sqrt(nu / x[0]) * ((x[1] * sin(v) * temp[0]) + (1 + x[1] * cos(v)) * temp[3]);//Vx
	result[4] = sqrt(nu / x[0]) * ((x[1] * sin(v) * temp[1]) + (1 + x[1] * cos(v)) * temp[4]);//Vy
	result[5] = sqrt(nu / x[0]) * temp(5);//Vz
	
	return 0;//exit
}

int main(int argc, char* argv[])
{
	X = (double *)malloc(N * sizeof(double));
	resultat = (double *)malloc(N * sizeof(double));
	X[0] =; 
	X[1] =;
	X[2] =;
	X[3] =;
	X[4] =;
	X[5] =;
	printmas(X, 6);
	kepler(x, resultat);
	printmas(resultat, 6);
	return 0;
}