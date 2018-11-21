#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS_GLOBALS
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <mpi.h>



// решение уравнения Ax=b методом Зейделя
// загрузка данных СЛУ из файла
int LoadMatrixFromFile(FILE *pFileIn, double *pCoefWork, double *pFreeMembWork, double *pSolWork, int inumRow, int inumCol);
// генерация псевдослучаных чисел в заданном интервале
double funRangedRand(int RangeMin, int RangeMax);
// расчет независисмой части метода Зейделя
double calcNDPart(double *pCoef, double *pSolx, int inumCol, int numstr);
// невязка = Ax - b
double funCalcDiscrepancy(double *pCoef, double *pFreeMemb, double *pSolx, int NumLine, int NumColum);




int main(int argc, char *argv[])
{
	int ProcNum; //кол-во процессов
	int ProcRank;//номер процесса
	MPI_Status Status;
	FILE *pFileIn;
	double *pMemCoef;                      // матрица A
	double *pMemFreeMemb;                  // вектор b
	double *pMemSolx;                      // вектор x             
	int  matrixsize;                       // кол-во элементов в матрице коэф-ов   
	int numCol;                        // количество столбцов (столбцы = строкам)
	int numRow;                        // количество строк
	double discr;                          // погрешность  
	double eps = 0.0000000001;
	int vret = 0;
	int RowiProc;                          // число обрабатываемых строк в процессе  
	double *pCoefVec;                      // передаваемая строка коффициентов
	double *pb;                            // передаваемый вектор b
	double *px;                            // передаваемый вектор x
	int blockelem;
	double xsi;
	double sum2;
	double sum1;
	double end = 1;
	double t1, t2, dt;
	double *pxw;                             // рабочий вектор



	MPI_Init(&argc, &argv);
	// определение количества процессов
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	// определение ранга процесса
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	pFileIn = fopen("input.txt", "r");
	if (pFileIn == NULL) { return -1; }

	// чтение размера матрицы (строки x столбцы)
	vret = fscanf(pFileIn, "%dx%d", &numRow, &numCol);
	if (vret <= 0) { return -2; }



	// выделение памяти для СЛУ
	pMemCoef = new double[numRow * numCol];
	pMemFreeMemb = new double[numRow];
	pMemSolx = new double[numRow];

	if (ProcRank == 0)
	{
		// загрузка СЛУ из файла
		LoadMatrixFromFile(pFileIn, pMemCoef, pMemFreeMemb, pMemSolx, numRow, numCol);
	}

	// определение количества строк, которое будет передаваться процессам
	// число строк кратно числу процессов
	RowiProc = numRow / ProcNum;
	blockelem = RowiProc * numCol;

	// для вектора с коэффициентами
	pCoefVec = new double[blockelem];
	// для вектора b
	pb = new double[RowiProc];
	// для вектора x
	px = new double[numRow];
	pxw = new double[numRow];//

	for (int vj = 0; vj < numRow; vj++)
	{
		pxw[vj] = 0;
		px[vj] = 0;
	}
	// рассылка всем процессам признака завершения работы программы
	MPI_Bcast(&end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// рассылка коэффициентов
	MPI_Scatter(pMemCoef, blockelem, MPI_DOUBLE, pCoefVec, blockelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// рассылка вектора b
	MPI_Scatter(pMemFreeMemb, RowiProc, MPI_DOUBLE, pb, RowiProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// рассылка вектора с приближенным решением 
	MPI_Bcast(pMemSolx, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pxw, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	dt = 0.0;

	// поиск решения происходит до тех пор, пока не будет достигнута заданная точность
	do
	{
		xsi = 0.0;
		t1 = MPI_Wtime();

		sum2 = 0.0;
		// расчет независимой части, вычисление суммы строк матрицы 
		sum2 = calcNDPart(pCoefVec, pMemSolx, numCol, ProcRank*RowiProc);
		sum1 = 0.0;
		// расчет зависимой части
		for (int vi = 0; vi < ProcRank; vi++)
		{
			for (int vj = 0; vj < RowiProc; vj++)
			{

			
				//printf("[%d] going to recv from %d proc with tag = %d\n", ProcRank, vi, vi*RowiProc + vj);
				//fflush(stdout);
				MPI_Recv(&xsi, 1, MPI_DOUBLE, vi, vi*RowiProc + vj, MPI_COMM_WORLD, &Status);
				pxw[vi*RowiProc + vj] = xsi;
				//printf("xsi%lf\n", xsi);
				// вычисление суммы строк матрицы L
				sum1 = sum1 + (pCoefVec[vi*RowiProc + vj] * pxw[vi*RowiProc + vj]);
				//printf("ProcRank=%d vi=%d vj=%d pCoefVec[vi*RowiProc + vj]=%lf  and  pxw[vi*RowiProc + vj]=%lf and sum1=%lf\n", ProcRank,vi, vj, pCoefVec[vi*RowiProc + vj], pxw[vi*RowiProc + vj], sum1);
				//fflush(stdout);
			}
		}
		// расчет приближенного решения
		pxw[ProcRank*RowiProc] = (pb[0] - sum1 - sum2) / pCoefVec[ProcRank*RowiProc];
		//printf("ProcRank=%d pb[0]=%lf pxw[ProcRank*RowiProc]=%lf sum1=%lf sum2=%lf  and  pCoefVec[ProcRank*RowiProc]=%lf and sum1=%lf\n", ProcRank, pb[0], pxw[ProcRank*RowiProc], sum1, sum2, pCoefVec[ProcRank*RowiProc], sum1);
		//fflush(stdout);

		// рассылка приближенного решения нижележащим процессам
		for (int vj = ProcRank + 1; vj < ProcNum; vj++)
		{
			//printf("Send-a:Proc=%d,Tag=%d,pxw[%d] = %lf\n", ProcRank, ProcRank*RowiProc, ProcRank*RowiProc, pxw[ProcRank*RowiProc]);
			//fflush(stdout);
			MPI_Send(&pxw[ProcRank*RowiProc], 1, MPI_DOUBLE, vj, ProcRank*RowiProc, MPI_COMM_WORLD);
		}

		for (int vk = 1; vk < RowiProc; vk++)
		{
			sum2 = 0.0;
			// расчет независимой части, вычисление суммы строк матрицы R
			sum2 = calcNDPart(pCoefVec + vk*numCol, pMemSolx, numCol, ProcRank*RowiProc + vk);
			sum1 = 0.0;
			// расчет зависимой части
			for (int vi = 0; vi < ProcRank*RowiProc + vk; vi++)
			{
				
				// на каждой итерации требуется предыдущее приближенное решение
				// прием приближенного рещения
				// вычисление суммы строк матрицы L
				sum1 = sum1 + (pCoefVec[vk*numCol + vi] * pxw[vi]);
			}

			// расчет приближенного решения
			pxw[ProcRank*RowiProc + vk] = (pb[vk] - sum1 - sum2) / pCoefVec[vk*numCol + vk + ProcRank*RowiProc];
			
			// рассылка приближенного решения нижележащим процессам
			for (int vj = ProcRank + 1; vj < ProcNum; vj++)
			{
				//printf("Send-b :Proc=%d,Tag=%d,pxw1[%d] = %lf\n", ProcRank, ProcRank*RowiProc + vk, ProcRank*RowiProc + vk, pxw[ProcRank*RowiProc + vk]);
				//fflush(stdout);
				MPI_Send(&pxw[ProcRank*RowiProc + vk], 1, MPI_DOUBLE, vj, ProcRank*RowiProc + vk, MPI_COMM_WORLD);
			}
		}

		// сбор приближенных решений на нулевом процессе
		MPI_Gather(&pxw[ProcRank*RowiProc], RowiProc, MPI_DOUBLE, px, RowiProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		if (ProcRank == 0)
		{

			// вычисление погрешности

			
			discr = fabs(px[0] - pMemSolx[0]);
			//printf("Discrepancy = %.12lf\n", discr);
			//fflush(stdout);
			// если достигнута заданная точность
			if (discr <= eps)
			{
				// завершение вычислений
				end = 0;
				// вывод результата
				/*for (int vs = 0; vs < numRow; vs++)
				{
					printf("x%d = %.12lf\n", vs, px[vs]);
				}*/
				/* невязка = Ax - b*/
				discr = funCalcDiscrepancy(pMemCoef, pMemFreeMemb, px, numRow, numCol);
				printf("Discrepancy = %.12lf\n", discr);

				printf("time = %lf\n", dt);

			}
			for (int vk = 0; vk < numRow; vk++)
			{
				pMemSolx[vk] = px[vk];
			}
		}
		// рассылка приближенного решения для дальнейших вычислений
		MPI_Bcast(pMemSolx, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(pxw, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		dt += (t2 - t1);
	} while (end);

	// освобождение памяти
	delete[] pMemCoef;
	delete[] pMemFreeMemb;
	delete[] pMemSolx;
	delete[] pCoefVec;
	delete[] pb;
	delete[] px;


	MPI_Finalize();


	return 0;
}



int LoadMatrixFromFile(FILE *pFileIn, double *pCoefWork, double *pFreeMembWork, double *pSolWork, int inumRow, int inumCol)
{

	int vret = 0;
	char vSym = '\0';


	// проверка количества сторок и столбцов
	if (inumRow <= 0 || inumCol <= 0 || inumRow != inumCol) { return -2; }
	// чтение матрицы  с коэф-ми
	for (int vi = 0; vi < inumRow; vi++)
	{
		vSym = '\0';
		for (int vj = 0; vj < inumCol; vj++)
		{
			vret = fscanf(pFileIn, "%lf", &pCoefWork[vi*inumRow + vj]);
			// если нет назначенных полей
			if (vret <= 0) { delete[]pCoefWork; pCoefWork = NULL; return -3; }
			fread(&vSym, sizeof(char), 1, pFileIn);
			// проверка формата матрицы коэф-ов
			if (vSym == ';')
			{
				// если количество столбцов не совпадает с реальным количеством столбцов
				if (inumCol - 1 != vj) { delete[] pCoefWork; pCoefWork = NULL; return -3; }
			}
		}
		// проверка формата матрицы коэф-ов
		if (vSym != ';') { delete[] pCoefWork; pCoefWork = NULL; return -3; }
	}
	// чтение свободных членов
	for (int vi = 0; vi < inumRow; vi++)
	{
		vSym = '\0';
		vret = fscanf(pFileIn, "%lf", &pFreeMembWork[vi]);
		// если нет назначенных полей
		if (vret <= 0) { delete[]pCoefWork; pCoefWork = NULL; delete[]pFreeMembWork; pFreeMembWork = NULL; return -4; }
		fread(&vSym, sizeof(char), 1, pFileIn);
		// если не соответствует формату матрицы свободных членов
		if (vSym != ';') { delete[]pCoefWork; pCoefWork = NULL; delete[]pFreeMembWork; pFreeMembWork = NULL; return -4; }
	}

	// заполнение вектора решения единицами
	for (int vi = 0; vi < inumRow; vi++)
	{
		pSolWork[vi] = 1.0;
	}

	fclose(pFileIn);


	return 0;
}

// расчет независисмой части метода Зейделя
double calcNDPart(double *pCoef, double *pSolx, int inumCol, int numstr)
{
	double vsum2 = 0.0;

	vsum2 = 0.0;
	for (int vk = numstr + 1; vk < inumCol; vk++)
	{
		vsum2 = vsum2 + pCoef[vk] * pSolx[vk];

		/*printf("vk=%d pCoef=%lf  and  pSolx=%lf and vsum2=%lf\n", vk,pCoef[vk], pSolx[vk], vsum2);
		fflush(stdout);*/
	}

	return vsum2;
}
// невязка = Ax - b
double funCalcDiscrepancy(double *pCoef, double *pFreeMemb, double *pSolx, int NumLine, int NumColum)
{
	int vi = 0;
	int vj = 0;
	double vsumAx = 0.0;            // сумма для матрицы Ax
	double vsumb = 0.0;             // сумма для матрицы b
	double vDesc = 0.0;

	for (vi = 0; vi < NumLine; vi++)
	{
		for (vj = 0; vj < NumColum; vj++)
		{
			// сумма для Ax
			vsumAx = vsumAx + (pCoef[vi*NumLine + vj] * pSolx[vj]);
		}
		// сумма для b
		vsumb = vsumb + pFreeMemb[vi];

		// деление полученных сумм на кол-во строк(столбцов)
		vDesc = fabs((vsumAx / NumLine) - (vsumb / NumLine));
	}

	return  vDesc;
}





