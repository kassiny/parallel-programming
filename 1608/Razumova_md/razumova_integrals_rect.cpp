#include <cmath>
#include <iostream>
#include <mpi.h>

using namespace std;

double func(double x, double y)
{
	return x * x + y * y;
}

double integral(double f(double, double), double x1, double y1, int nx, int ny, double hx, double hy)
{
	double sum = 0;
	double x, y;

	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			x = x1 + i * hx + hx / 2.0;
			y = y1 + j * hy + hy / 2.0;

			sum += f(x, y)*hx*hy;
		}
	}
	cout << "x1 " << x1 << " y1 " << y1 << " nx " << nx << " hx " << hx << " hy " << hy << " mysum " << sum << endl;
	return sum;
}


#define LB_X    0
#define UB_X    10
#define LB_Y    0
#define UB_Y    10
#define N_TOTAL 1000.0


int main(int argc, char *argv[])
{
	int myid, np, n_per_p;
	double sum, total;
	double hx, hy, my_range;
	double start_x, start_y;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	hx = float((UB_X - LB_X) / N_TOTAL);
	hy = float((UB_Y - LB_Y) / N_TOTAL);
	start_x = LB_X + (double)myid*((UB_X - LB_X) / (double)np);

	n_per_p = N_TOTAL / np;
	if (myid == np - 1) {
		n_per_p += N_TOTAL - n_per_p * np;
	}

	sum = integral(func, start_x, LB_Y, n_per_p, N_TOTAL, hx, hy);
	//cout << myid << ": Start " << start_x << " " <<  start_y << "NUM " << N_TOTAL/np << endl;

	if (myid == 0) {
		MPI_Reduce(&sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		cout << "Got " << total << endl;
	}
	else {
		MPI_Reduce(&sum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	cout << integral(func,LB_X, LB_Y, N_TOTAL, N_TOTAL, float((UB_X - LB_X) / N_TOTAL), float((UB_Y - LB_Y) / N_TOTAL)) << endl;
	return 0;
}
