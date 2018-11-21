#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS_GLOBALS
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <mpi.h>



// ������� ��������� Ax=b ������� �������
// �������� ������ ��� �� �����
int LoadMatrixFromFile(FILE *pFileIn, double *pCoefWork, double *pFreeMembWork, double *pSolWork, int inumRow, int inumCol);
// ��������� �������������� ����� � �������� ���������
double funRangedRand(int RangeMin, int RangeMax);
// ������ ������������ ����� ������ �������
double calcNDPart(double *pCoef, double *pSolx, int inumCol, int numstr);
// ������� = Ax - b
double funCalcDiscrepancy(double *pCoef, double *pFreeMemb, double *pSolx, int NumLine, int NumColum);




int main(int argc, char *argv[])
{
	int ProcNum; //���-�� ���������
	int ProcRank;//����� ��������
	MPI_Status Status;
	FILE *pFileIn;
	double *pMemCoef;                      // ������� A
	double *pMemFreeMemb;                  // ������ b
	double *pMemSolx;                      // ������ x             
	int  matrixsize;                       // ���-�� ��������� � ������� ����-��   
	int numCol;                        // ���������� �������� (������� = �������)
	int numRow;                        // ���������� �����
	double discr;                          // �����������  
	double eps = 0.0000000001;
	int vret = 0;
	int RowiProc;                          // ����� �������������� ����� � ��������  
	double *pCoefVec;                      // ������������ ������ ������������
	double *pb;                            // ������������ ������ b
	double *px;                            // ������������ ������ x
	int blockelem;
	double xsi;
	double sum2;
	double sum1;
	double end = 1;
	double t1, t2, dt;
	double *pxw;                             // ������� ������



	MPI_Init(&argc, &argv);
	// ����������� ���������� ���������
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	// ����������� ����� ��������
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	pFileIn = fopen("input.txt", "r");
	if (pFileIn == NULL) { return -1; }

	// ������ ������� ������� (������ x �������)
	vret = fscanf(pFileIn, "%dx%d", &numRow, &numCol);
	if (vret <= 0) { return -2; }



	// ��������� ������ ��� ���
	pMemCoef = new double[numRow * numCol];
	pMemFreeMemb = new double[numRow];
	pMemSolx = new double[numRow];

	if (ProcRank == 0)
	{
		// �������� ��� �� �����
		LoadMatrixFromFile(pFileIn, pMemCoef, pMemFreeMemb, pMemSolx, numRow, numCol);
	}

	// ����������� ���������� �����, ������� ����� ������������ ���������
	// ����� ����� ������ ����� ���������
	RowiProc = numRow / ProcNum;
	blockelem = RowiProc * numCol;

	// ��� ������� � ��������������
	pCoefVec = new double[blockelem];
	// ��� ������� b
	pb = new double[RowiProc];
	// ��� ������� x
	px = new double[numRow];
	pxw = new double[numRow];//

	for (int vj = 0; vj < numRow; vj++)
	{
		pxw[vj] = 0;
		px[vj] = 0;
	}
	// �������� ���� ��������� �������� ���������� ������ ���������
	MPI_Bcast(&end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// �������� �������������
	MPI_Scatter(pMemCoef, blockelem, MPI_DOUBLE, pCoefVec, blockelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// �������� ������� b
	MPI_Scatter(pMemFreeMemb, RowiProc, MPI_DOUBLE, pb, RowiProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// �������� ������� � ������������ �������� 
	MPI_Bcast(pMemSolx, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pxw, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	dt = 0.0;

	// ����� ������� ���������� �� ��� ���, ���� �� ����� ���������� �������� ��������
	do
	{
		xsi = 0.0;
		t1 = MPI_Wtime();

		sum2 = 0.0;
		// ������ ����������� �����, ���������� ����� ����� ������� 
		sum2 = calcNDPart(pCoefVec, pMemSolx, numCol, ProcRank*RowiProc);
		sum1 = 0.0;
		// ������ ��������� �����
		for (int vi = 0; vi < ProcRank; vi++)
		{
			for (int vj = 0; vj < RowiProc; vj++)
			{

			
				//printf("[%d] going to recv from %d proc with tag = %d\n", ProcRank, vi, vi*RowiProc + vj);
				//fflush(stdout);
				MPI_Recv(&xsi, 1, MPI_DOUBLE, vi, vi*RowiProc + vj, MPI_COMM_WORLD, &Status);
				pxw[vi*RowiProc + vj] = xsi;
				//printf("xsi%lf\n", xsi);
				// ���������� ����� ����� ������� L
				sum1 = sum1 + (pCoefVec[vi*RowiProc + vj] * pxw[vi*RowiProc + vj]);
				//printf("ProcRank=%d vi=%d vj=%d pCoefVec[vi*RowiProc + vj]=%lf  and  pxw[vi*RowiProc + vj]=%lf and sum1=%lf\n", ProcRank,vi, vj, pCoefVec[vi*RowiProc + vj], pxw[vi*RowiProc + vj], sum1);
				//fflush(stdout);
			}
		}
		// ������ ������������� �������
		pxw[ProcRank*RowiProc] = (pb[0] - sum1 - sum2) / pCoefVec[ProcRank*RowiProc];
		//printf("ProcRank=%d pb[0]=%lf pxw[ProcRank*RowiProc]=%lf sum1=%lf sum2=%lf  and  pCoefVec[ProcRank*RowiProc]=%lf and sum1=%lf\n", ProcRank, pb[0], pxw[ProcRank*RowiProc], sum1, sum2, pCoefVec[ProcRank*RowiProc], sum1);
		//fflush(stdout);

		// �������� ������������� ������� ����������� ���������
		for (int vj = ProcRank + 1; vj < ProcNum; vj++)
		{
			//printf("Send-a:Proc=%d,Tag=%d,pxw[%d] = %lf\n", ProcRank, ProcRank*RowiProc, ProcRank*RowiProc, pxw[ProcRank*RowiProc]);
			//fflush(stdout);
			MPI_Send(&pxw[ProcRank*RowiProc], 1, MPI_DOUBLE, vj, ProcRank*RowiProc, MPI_COMM_WORLD);
		}

		for (int vk = 1; vk < RowiProc; vk++)
		{
			sum2 = 0.0;
			// ������ ����������� �����, ���������� ����� ����� ������� R
			sum2 = calcNDPart(pCoefVec + vk*numCol, pMemSolx, numCol, ProcRank*RowiProc + vk);
			sum1 = 0.0;
			// ������ ��������� �����
			for (int vi = 0; vi < ProcRank*RowiProc + vk; vi++)
			{
				
				// �� ������ �������� ��������� ���������� ������������ �������
				// ����� ������������� �������
				// ���������� ����� ����� ������� L
				sum1 = sum1 + (pCoefVec[vk*numCol + vi] * pxw[vi]);
			}

			// ������ ������������� �������
			pxw[ProcRank*RowiProc + vk] = (pb[vk] - sum1 - sum2) / pCoefVec[vk*numCol + vk + ProcRank*RowiProc];
			
			// �������� ������������� ������� ����������� ���������
			for (int vj = ProcRank + 1; vj < ProcNum; vj++)
			{
				//printf("Send-b :Proc=%d,Tag=%d,pxw1[%d] = %lf\n", ProcRank, ProcRank*RowiProc + vk, ProcRank*RowiProc + vk, pxw[ProcRank*RowiProc + vk]);
				//fflush(stdout);
				MPI_Send(&pxw[ProcRank*RowiProc + vk], 1, MPI_DOUBLE, vj, ProcRank*RowiProc + vk, MPI_COMM_WORLD);
			}
		}

		// ���� ������������ ������� �� ������� ��������
		MPI_Gather(&pxw[ProcRank*RowiProc], RowiProc, MPI_DOUBLE, px, RowiProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		if (ProcRank == 0)
		{

			// ���������� �����������

			
			discr = fabs(px[0] - pMemSolx[0]);
			//printf("Discrepancy = %.12lf\n", discr);
			//fflush(stdout);
			// ���� ���������� �������� ��������
			if (discr <= eps)
			{
				// ���������� ����������
				end = 0;
				// ����� ����������
				/*for (int vs = 0; vs < numRow; vs++)
				{
					printf("x%d = %.12lf\n", vs, px[vs]);
				}*/
				/* ������� = Ax - b*/
				discr = funCalcDiscrepancy(pMemCoef, pMemFreeMemb, px, numRow, numCol);
				printf("Discrepancy = %.12lf\n", discr);

				printf("time = %lf\n", dt);

			}
			for (int vk = 0; vk < numRow; vk++)
			{
				pMemSolx[vk] = px[vk];
			}
		}
		// �������� ������������� ������� ��� ���������� ����������
		MPI_Bcast(pMemSolx, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(pxw, numRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		dt += (t2 - t1);
	} while (end);

	// ������������ ������
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


	// �������� ���������� ������ � ��������
	if (inumRow <= 0 || inumCol <= 0 || inumRow != inumCol) { return -2; }
	// ������ �������  � ����-��
	for (int vi = 0; vi < inumRow; vi++)
	{
		vSym = '\0';
		for (int vj = 0; vj < inumCol; vj++)
		{
			vret = fscanf(pFileIn, "%lf", &pCoefWork[vi*inumRow + vj]);
			// ���� ��� ����������� �����
			if (vret <= 0) { delete[]pCoefWork; pCoefWork = NULL; return -3; }
			fread(&vSym, sizeof(char), 1, pFileIn);
			// �������� ������� ������� ����-��
			if (vSym == ';')
			{
				// ���� ���������� �������� �� ��������� � �������� ����������� ��������
				if (inumCol - 1 != vj) { delete[] pCoefWork; pCoefWork = NULL; return -3; }
			}
		}
		// �������� ������� ������� ����-��
		if (vSym != ';') { delete[] pCoefWork; pCoefWork = NULL; return -3; }
	}
	// ������ ��������� ������
	for (int vi = 0; vi < inumRow; vi++)
	{
		vSym = '\0';
		vret = fscanf(pFileIn, "%lf", &pFreeMembWork[vi]);
		// ���� ��� ����������� �����
		if (vret <= 0) { delete[]pCoefWork; pCoefWork = NULL; delete[]pFreeMembWork; pFreeMembWork = NULL; return -4; }
		fread(&vSym, sizeof(char), 1, pFileIn);
		// ���� �� ������������� ������� ������� ��������� ������
		if (vSym != ';') { delete[]pCoefWork; pCoefWork = NULL; delete[]pFreeMembWork; pFreeMembWork = NULL; return -4; }
	}

	// ���������� ������� ������� ���������
	for (int vi = 0; vi < inumRow; vi++)
	{
		pSolWork[vi] = 1.0;
	}

	fclose(pFileIn);


	return 0;
}

// ������ ������������ ����� ������ �������
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
// ������� = Ax - b
double funCalcDiscrepancy(double *pCoef, double *pFreeMemb, double *pSolx, int NumLine, int NumColum)
{
	int vi = 0;
	int vj = 0;
	double vsumAx = 0.0;            // ����� ��� ������� Ax
	double vsumb = 0.0;             // ����� ��� ������� b
	double vDesc = 0.0;

	for (vi = 0; vi < NumLine; vi++)
	{
		for (vj = 0; vj < NumColum; vj++)
		{
			// ����� ��� Ax
			vsumAx = vsumAx + (pCoef[vi*NumLine + vj] * pSolx[vj]);
		}
		// ����� ��� b
		vsumb = vsumb + pFreeMemb[vi];

		// ������� ���������� ���� �� ���-�� �����(��������)
		vDesc = fabs((vsumAx / NumLine) - (vsumb / NumLine));
	}

	return  vDesc;
}





