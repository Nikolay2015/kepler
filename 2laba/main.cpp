#include <iostream>//#include <stdio.h>
#include <math.h>
#include <stdio.h>
#define M_PI	3.14159265358979323846

unsigned int n = 6;// int n = 6;
double *r;//x, y, z, vx, vy, vz
double *result;

using namespace std;

double V_(double c[3], double f[3], double *x, double f_, double r_)
{
	if(((f[1] * x[2] - f[2] * x[1]) * c[0] + (f[2] * x[0] - f[0] * x[2]) * c[1] + (f[0] * x[1] - f[1] * x[0]) * c[2]) > 0)
		return (acos((f[0] * x[0] + f[1] * x[1] + f[2] * x[2]) / (f_ * r_)));
	else
		return (2 * M_PI - acos((f[0] * x[0] + f[1] * x[1] + f[2] * x[2]) / (f_ * r_)));
}

void obrkepler(double *X, double t)
{
	double R;
	double V;
	double H;
	//double mu = 398600.0;//km^3/sec;
	double mu = 398600e9;
	double C[3];// cx, cy, cz
	double F[3];// fx, fy, fz
	double F_;// sqrt(fx+fy+fz)
	double C_;//sqrt(cx+cy+cz)
	double C_2;//sqrt(cx + cy)
	double E;//
	double a;//p/1-e2
	double n;//sqrt(mu / a^3)
	
	R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	V = sqrt(X[3] * X[3] + X[4] * X[4] + X[5] * X[5]);
	H = V * V - ((2.0 * mu)/R);
	C[0] = X[1] * X[5] - X[2] * X[4];
	C[1] = X[2] * X[3] - X[0] * X[5];
	C[2] = X[0] * X[4] - X[1] * X[3];
	C_ = sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
	result[0] = (C_ * C_) / mu;//p
	F[0] = X[4] * C[2] - X[5] * C[1] - (mu / R) * X[0];//
	F[1] = X[5] * C[0] - X[3] * C[3] - (mu / R) * X[1];
	F[2] = X[3] * C[1] - X[4] * C[0] - (mu / R) * X[2];
	F_ = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
	result[1] = F_ / mu;// e
	result[4] = acos(C[2] / C_);// i
	C_2 = sqrt(C[0] * C[0] + C[1] * C[1]);
	if ((C[0] / C_2) > 0)
		result[5] = acos(-(C[1] / C_2));//O
	else
		result[5] = 2.0 * M_PI - acos(-(C[1] / C_2));//O
	double temp;
	temp = (1.0 / C_2) 
	* (F[2] * (C[0] * C[0] + C[1] * C[1]) - C[2] * (F[0] * C[0] + F[1] * C[1]));
	
	if (temp  > 0)
		result[2] = acos(((-C[1] / C_2) * F[0] + (C[0] / C_2) * F[1]) / F_);// omeg
	else
		result[2] = 2 * M_PI 
	- acos(((-C[1] / C_2) * F[0] + (C[0] / C_2) * F[1]) / F_);// omeg
	
	if (H < 0)//eleptic orbit
	{
		E = 2 * atan(sqrt((1.0 - result[1]) / (1.0 + result[1])) * tan(V_(C, F, X, F_, C_) / 2.0));
		a = result[0] / (1 - result[1] * result[1]);
		n = sqrt(mu / (a * a * a));
		result[3] = t - ((E - (result[1] * sin(E))) / n);
		//result[3] = n * (t - result[3]);
	}
	else if( H = 0)//parabolitic orbit
		{
			temp = tan(V_(C, F, X, F_, C_) / 2.0) + (1.0 / 3.0) * pow(tan(V_(C, F, X, F_, C_) / 2.0), 3);
			result[3] = t - ((1.0 / 2.0) * sqrt(pow(result[0], 3) / mu)) * temp;
			//result[3] = n * (t - result[3]);
		}
		else//giperolitic orbit
		{
			temp = 2.0 * atanh(sqrt((1 - result[1]) / (1 + result[1])) * tan(V_(C, F, X, F_, C_) / 2.0));
			result[3] = t - (sqrt(pow(result[0] / (result[1] * result[1] - 1.0), 3) / mu) * (result[1] * sin(temp) - temp));
			//result[3] = n * (t - result[3]);
		}
	
}

void print(double *x, unsigned int n)
{
	for(unsigned int i = 0; i < n; i++)//unsigned int i; for (i = 0; i < n; i++)
	{
		printf("%f\n", x[i]);
	}
}
void print(double *x, unsigned int n, int scilab)
{
	printf("a %f\n", x[0]);
	printf("e %f\n", x[1]);
	printf("i %f\n", x[4]);
	printf("omega %f\n", x[2]);
	printf("O %f\n", x[5]);
	printf("M %f\n", x[3]);
}
int main()
{
	r = new double[n];// r = (double *)malloc(n * sizeof(double));
	result = new double[n];//p, e, omeg, tau, i, O
	r[0] = 7000e3;
	r[1] = 1000e3;
	r[2] = -500e3;
	r[3] = 1e3;
	r[4] = 2e3;
	r[5] = 7e3;
	print(r, 6);
	obrkepler(r, 0);
	print(result, 6, 1);
	//r[0] = //
	return 0;
}
