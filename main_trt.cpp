#pragma warning(disable:4996)
#include <iostream>
#include<mpi.h>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <iomanip>
#include <fstream>
#include<cstring>
#include <sstream>
#include"MPI_Def.h"

using namespace std;

#define i_min 1
#define j_min 1
#define TIM   4

// MPI data structure 
MPI_Comm comm1;
MPI_Status Status1;
MPI_Datatype SendRecvDerive;

// coordinate data structure 
int grid[CRAT_DIM];         // x,y方向的网格尺度
int MatScale[CRAT_DIM];     // 每一小块的矩阵的大小
int coords[CRAT_DIM];       // 每一个rank的坐标
int Dims[CRAT_DIM];         // 每一个维度的块的个数
int periods[CRAT_DIM];      // 笛卡尔坐标是否具有周期性
int Disp[CRAT_DIM];         // 每一个子块首个元素的全局坐标
int Neighbor[4];            // 左右，上下邻居
int i_max, j_max;
int i_min_size, j_min_size, i_max_size, j_max_size;

//this struct is used to send or receive the data
struct MPI_SR_Data
{
	double F_send[9], u_send[2], rho_send;
};
// compute data structure 
double         ***ff, ***F, ***f, ***u, ***u0, **rho, ***Drag_F, Force;
double        fx = 0.0, fy = 0.0, temp_1=0.0, temp_2=0.0,totaltem1=0.0,totaltem2=0.0;
const double  c = 1.0;
double        flag = 0.0;
double        totalflag = 0.0;
int           **Soid;
MPI_SR_Data   *FDownRecv, *FUpRecv, *FLeftRecv, *FRightRecv,
              *FDownSend, *FUpSend, *FLeftSend, *FRightSend;

double       Re, dx, dy, Lx, Ly, dt, error = 2.0, Sniu, S_q;
const double rho0 = 1.0;
int        NX, NY, ip = 4, jp = 2, i, j;
double theta = 0.0*PI;//3.1415926 / 5.0;
double     Volume_Fraction;
double     totalerror = 0.0, totalfx = 0.0, totalfy = 0.0;
double         w[Q] = { 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 };
const int e_f[Q][2] = { { 0,0 },{ 1,0 },{ 0,1 },{ -1,0 },{ 0,-1 },{ 1,1 },{ -1,1 },{ -1,-1 },{ 1,-1 } };    // 9 direction
const int ne_f[Q] = { 0,3,4,1,2,7,8,5,6 };
int   break_while = 1;

void     Space_Malloc(int BlockSize[2]);
void     Space_Free(int BlockSize[2]);
void     Set_Coord_Info();
void     Init_Field();
void     Soid_Flag();
double   feq(int k, double rho, double u[2]);
void     Commucate_Data();
void     Evolution();
void     Error();
void     Output_Result(int myrank);
void    Drag_Force();


int main(int argc, char *argv[])

{
	
	MPI_Init(&argc, &argv);
	int myrank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	Set_Coord_Info();
	//cout << "this rank :" << myrank << " does not aborted" << endl;
	Init_Field();
	Error();
	int n = 0;
	/*cout << "this rank :" << myrank << " does not aborted" << endl;*/
	while (break_while)//break_while
	{
		Evolution();


		if (n % 100 == 0)
		{

			Error();
			totalerror = 0.0;
			//cout << myrank<<" th the error is :" << error << endl;
		
			MPI_Reduce(&temp_1, &totaltem1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&temp_2, &totaltem2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			totalerror = totaltem1/(totaltem2+1.0e-30);
		}
		if (myrank == 0 && n % 100 == 0)
		{
			printf("the computation step is : %.8d \n", n);

			printf("the error of total error  is : %.12lf !!\n", totalerror);
			//cout << "the speed is :" << u[1][(j_min + j_max_size) / 2][0]- u0[1][(j_min + j_max_size) / 2][0] << endl;
		}

		if (n % 1000 == 0)
		{
			MPI_Bcast(&totalerror, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			MPI_Barrier(MPI_COMM_WORLD);
			if (totalerror<1.0e-11)
			{
				break_while = 0;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		n++;
		
	}	
	Drag_Force();
    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&fx, &totalfx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&fy, &totalfy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (myrank == 0)
	{
		printf("the force of total force in x is : %.12lf !!\n", totalfx);
		printf("the force of total force in y is : %.12lf !!\n", totalfy);
	}
	//printf("the rank %2d force of total force in x is : %.12lf !!\n",myrank, totalfx);
	Output_Result(myrank);
	Space_Free(MatScale);


	MPI_Finalize();
	system("pause");
	return 0;
}

void    Space_Malloc(int BlockSize[2])
{
	
	i_min_size = 1, i_max_size = BlockSize[0] + 1;
	j_min_size = 1, j_max_size = BlockSize[1] + 1;

	ALLOC_2d_array(double, rho, BlockSize[0] + 2, BlockSize[1] + 2);
	ALLOC_2d_array(int, Soid, BlockSize[0] + 2, BlockSize[1] + 2);
	ALLOC_3d_array(double, ff, BlockSize[0] + 2, BlockSize[1] + 2, Q);
	ALLOC_3d_array(double, f, BlockSize[0] + 2, BlockSize[1] + 2, Q);
	ALLOC_3d_array(double, F, BlockSize[0] + 2, BlockSize[1] + 2, Q);
	ALLOC_3d_array(double, u, BlockSize[0] + 2, BlockSize[1] + 2, CRAT_DIM);
	ALLOC_3d_array(double, u0, BlockSize[0] + 2, BlockSize[1] + 2, CRAT_DIM);
	ALLOC_3d_array(double, Drag_F, BlockSize[0] + 2, BlockSize[1] + 2, CRAT_DIM);
	FLeftSend = new MPI_SR_Data[BlockSize[1]], FRightSend = new MPI_SR_Data[BlockSize[1]];
	FLeftRecv = new MPI_SR_Data[BlockSize[1]], FRightRecv = new MPI_SR_Data[BlockSize[1]];

	FUpSend = new MPI_SR_Data[BlockSize[0] + 2], FDownSend = new MPI_SR_Data[BlockSize[0] + 2];
	FUpRecv = new MPI_SR_Data[BlockSize[0] + 2], FDownRecv = new MPI_SR_Data[BlockSize[0] + 2];
}



/************************************************************************/
/* this is used to free the space                                       */
/************************************************************************/
void    Space_Free(int BlockSize[2])
{

	FREE_2d_array(double, rho, BlockSize[0] + 2);
	FREE_2d_array(int, Soid, BlockSize[0] + 2);
	FREE_3d_array(double, f, BlockSize[0] + 2, BlockSize[1] + 2);
	FREE_3d_array(double, ff, BlockSize[0] + 2, BlockSize[1] + 2);
	FREE_3d_array(double, F, BlockSize[0] + 2, BlockSize[1] + 2);
	FREE_3d_array(double, u, BlockSize[0] + 2, BlockSize[1] + 2);
	FREE_3d_array(double, u0, BlockSize[0] + 2, BlockSize[1] + 2);
	FREE_3d_array(double, Drag_F, BlockSize[0] + 2, BlockSize[1] + 2);

	delete[]FLeftSend, delete[]FRightSend;
	delete[]FLeftRecv, delete[]FRightRecv;
	delete[]FUpSend, delete[]FDownSend;
	delete[]FUpRecv, delete[]FDownRecv;

}

void    Set_Coord_Info()
{
	int rank, procs;
	coords[0] = 0;
	coords[1] = 0;
	Dims[0] = Dims[1] = 0;

	// set the period boundary

	periods[0] = periods[1] = 1;


	//get the number of the procs and the rank number
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	if (rank == 0)
	{
		cout << "please input the Grid points in x direction:" << endl;
		cin >> grid[0];
		cout << "please input the Grid points in y direction : " << endl;
		cin >> grid[1];
		cout << "please input the volume fraction of the flow: " << endl;
		cin >> Volume_Fraction;

	}



	MPI_Bcast(grid, CRAT_DIM, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Volume_Fraction, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	theta = Volume_Fraction*PI;
	NX = grid[0], NY = grid[1];
	MPI_Dims_create(procs, CRAT_DIM, Dims);
	MPI_Cart_create(MPI_COMM_WORLD, CRAT_DIM, Dims, periods, 1, &comm1);

	//get the coordinate from the rank 
	MPI_Cart_coords(comm1, rank, CRAT_DIM, coords);

	//每个子块的大小
	MatScale[0] = BLOCK_SIZE(coords[0], Dims[0], grid[0]);
	MatScale[1] = BLOCK_SIZE(coords[1], Dims[1], grid[1]);


	i_max_size = MatScale[0] + 1;
	j_max_size = MatScale[1] + 1;
	//计算每一个子块的偏移量
	Disp[0] = BLOCK_LOW(coords[0], Dims[0], grid[0]);
	Disp[1] = BLOCK_LOW(coords[1], Dims[1], grid[1]);

	Space_Malloc(MatScale);


	int RightCoords[2], LeftCoords[2], UpperCoords[2], DownCoords[2];
	LeftCoords[0] = coords[0] - 1;     LeftCoords[1] = coords[1];
	RightCoords[0] = coords[0] + 1;     RightCoords[1] = coords[1];
	UpperCoords[0] = coords[0];     UpperCoords[1] = coords[1] + 1;
	DownCoords[0] = coords[0];     DownCoords[1] = coords[1] - 1;



	MPI_Cart_rank(comm1, LeftCoords, Neighbor);      //左边进程的编号值
	MPI_Cart_rank(comm1, RightCoords, Neighbor + 1);    //右边进程的编号值
	MPI_Cart_rank(comm1, UpperCoords, Neighbor + 2);  //上边进程的编号值
	MPI_Cart_rank(comm1, DownCoords, Neighbor + 3);  //下边进程的编号值

}
void    Init_Field()
{
	int myrank;
	dx = 1.0;// double(NX);
	dy = 1.0;// double(NX);
	Lx = dx*double(NX-2);
	Ly = dy*double(NY-2);



	dt = dx;

	Force = 0.5 / Lx / Ly;
	double niu = 1.0;
	Sniu = -2.0* dt / (6.0*niu + dt);
	S_q  = -8.0*(2 + Sniu) / (8 + Sniu); //!!!!!
	Soid_Flag();
	for (int ii = i_min; ii < i_max_size; ii++)
	{
		for (int jj =j_min; jj < j_max_size; jj++)
		{
			u[ii][jj][0] = 0.00;
			u[ii][jj][1] = 0.00;
			rho[ii][jj] =  rho0;
			if (Soid[ii][jj] == 1)
				rho[ii][jj] = 0.0;
			for (int k = 0; k < Q; k++)
			{
				f [ii][jj][k] = feq(k, rho[ii][jj], u[ii][jj]);
			}
		}
	}
}
double  feq(int k1, double rho1, double u[2])   //  equilibrium distribution
{
	double eu, uv, feq;
	eu = (e_f[k1][0] * u[0] + e_f[k1][1] * u[1]);
	uv = (u[0] * u[0] + u[1] * u[1]);
	feq = w[k1] * rho1*(1.0 + 3.0*eu / c + 4.5*eu*eu / (c*c) - 1.5*uv / (c*c));
	return feq;
}
/*void    Soid_Flag()
{
	int *flag_left_send, *flag_right_send, *flag_up_send, *flag_down_send;
	int *flag_left_recv, *flag_right_recv, *flag_up_recv, *flag_down_recv;
	// initial the array to be sent
	flag_left_send  = new int[MatScale[1]];  flag_up_send   = new int[MatScale[0] + CRAT_DIM];
	flag_right_send = new int[MatScale[1]];  flag_down_send = new int[MatScale[0] + CRAT_DIM];
	flag_left_recv  = new int[MatScale[1]];  flag_up_recv   = new int[MatScale[0] + CRAT_DIM];
	flag_right_recv = new int[MatScale[1]];  flag_down_recv = new int[MatScale[0] + CRAT_DIM];


	int i_tid, j_tid;
	int i_dis, j_dis;
	i_tid = Disp[0], j_tid = Disp[1];
	for (int i = i_min; i<i_max_size; i++)
	{
		i_tid = Disp[0] + i - 1;
		for (int j = j_min; j<j_max_size; j++)
		{
			Soid[i][j] = 0;
			j_tid = Disp[1] + j - 1;
			if (i_tid + 1>NX / 2)
			{
				i_dis = (NX - i_tid)*TIM;
			}
			else
				i_dis = fabs(i_tid + 1)*TIM;
			if (j_tid + 1>NY / 2)
			{
				j_dis = (NY - j_tid)*TIM;
			}
			else
				j_dis = fabs(j_tid + 1)*TIM;
			if (i_dis <= NX&&j_dis <= NY)
			{
				Soid[i][j] = 1;
				if (i_dis == NX || j_dis == NY)
				{
					Soid[i][j] = 2;
				}
			}

		}
	}



	//填充左右数据
	for (int i = 0; i < MatScale[1]; i++)
	{
		flag_left_send [i] = Soid[i_min        ][i + 1];
		flag_right_send[i] = Soid[i_max_size -1][i + 1];
	}
	//数据从左往右传递，
	MPI_Sendrecv(flag_right_send, MatScale[1],   MPI_INT, Neighbor[1], 1,
			     flag_left_recv , MatScale[1],   MPI_INT, Neighbor[0], 1,MPI_COMM_WORLD,&Status1);



	//数据从右往左传递，
	MPI_Sendrecv(flag_left_send , MatScale[1],   MPI_INT, Neighbor[0], 1,
		         flag_right_recv, MatScale[1],   MPI_INT, Neighbor[1], 1,MPI_COMM_WORLD,&Status1);

	//将传递的数据进行解包
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < MatScale[1]; i++)
	{
		Soid[0         ][i + 1]= flag_left_recv [i];
		Soid[i_max_size][i + 1]= flag_right_recv[i];
	}

	
	//填充上下数据
	for (int j = 0; j < MatScale[0]+2; j++)
	{
		flag_up_send  [j] = Soid[j][j_max_size-1];
		flag_down_send[j] = Soid[j][j_min       ];
	}

	//数据从上往下传递，
	MPI_Sendrecv(flag_down_send, MatScale[0] + 2, MPI_INT, Neighbor[3], 1,
		         flag_up_recv  , MatScale[0] + 2, MPI_INT, Neighbor[2], 1, MPI_COMM_WORLD, &Status1);
	//数据从下往上传递

	MPI_Sendrecv(flag_up_send  , MatScale[0] + 2, MPI_INT, Neighbor[2], 1,
		         flag_down_recv, MatScale[0] + 2, MPI_INT, Neighbor[3], 1, MPI_COMM_WORLD, &Status1);

	//填充上下的数据
	for (int j = 0; j < MatScale[0] + 2; j++)
	{
		Soid[j][j_max_size] = flag_up_recv  [j];
		Soid[j][0]          = flag_down_recv[j];
	}


	delete []flag_left_send  ;      delete[]flag_up_send   ;
	delete []flag_right_send ;      delete[]flag_down_send ;
	delete []flag_left_recv  ;      delete[]flag_up_recv   ;
	delete []flag_right_recv ;      delete[]flag_down_recv ;

}*/

//void    Commucate_Data()
//{
//	//填充左右数据
//	int k = 0;
//	for (int i = 0; i < MatScale[1]; i++)
//	{
//		FLeftSend [i]. rho_send  = rho[i_min][i + 1];
//		FLeftSend [i]. u_send[0] = u  [i_min][i + 1][0];
//		FLeftSend [i]. u_send[1] = u  [i_min][i + 1][1];
//
//		FRightSend[i].rho_send  = rho[i_max_size - 1][i + 1];
//		FRightSend[i].u_send[0] = u  [i_max_size - 1][i + 1][0];
//		FRightSend[i].u_send[1] = u  [i_max_size - 1][i + 1][1];
//		for ( k = 0; k < Q; k++)
//		{
//			FRightSend[i].F_send [k] = F [i_max_size - 1][i + 1][k];
//			FRightSend[i].f_send [k] = f [i_max_size - 1][i + 1][k];
//			FRightSend[i].ff_send[k] = ff[i_max_size - 1][i + 1][k];
//			FLeftSend [i].F_send [k] = F [i_min][i + 1][k];
//			FLeftSend [i].f_send [k] = f [i_min][i + 1][k];
//			FLeftSend [i].ff_send[k] = ff[i_min][i + 1][k];
//		}
//		
//	}
//	//数据从左往右传递，
//	MPI_Send(FLeftSend /*UleftSend*/, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[0], 99, comm1);
//	MPI_Recv(FRightRecv/*UrightRec*/, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[1], 99, comm1, &Status1);
//
//	MPI_Send(FRightSend/*UrightSend*/, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[1], 99, comm1);
//	MPI_Recv(FLeftRecv /*UleftRec*/,   MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[0], 99, comm1, &Status1);
//	//MPI_Sendrecv(FRightSend, MatScale[1]*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[1], 1,
//	//	         FLeftRecv , MatScale[1]*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[0], 1, MPI_COMM_WORLD, &Status1);
//
//
//
//	////数据从右往左传递，
//	//MPI_Sendrecv(FLeftSend , MatScale[1]*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[0], 1,
//	//	         FRightRecv, MatScale[1]*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[1], 1, MPI_COMM_WORLD, &Status1);
//
//	//将传递的数据进行解包
//	MPI_Barrier(MPI_COMM_WORLD);
//	for (int i = 0; i < MatScale[1]; i++)
//	{
//		rho[0][i + 1]          = FLeftRecv[i].rho_send;
//		u  [0][i + 1][0]       = FLeftRecv[i].u_send[0];
//		u  [0][i + 1][1]       = FLeftRecv[i].u_send[1];
//
//		rho[i_max_size][i + 1]    = FRightRecv[i].rho_send;
//		u  [i_max_size][i + 1][0] = FRightRecv[i].u_send[0];
//		u  [i_max_size][i + 1][1] = FRightRecv[i].u_send[1];
//		for ( k = 0; k < Q; k++)
//		{
//			F [0][i + 1][k]          = FLeftRecv [i]. F_send[k];
//			//f [0][i + 1][k]          = FLeftRecv [i].f_send [k];
//			//ff[0][i + 1][k]          = FLeftRecv [i].ff_send[k];
//			F[i_max_size][i + 1][k] = FRightRecv[i].F_send[k];
//			//f[i_max_size][i + 1][k] = FRightRecv[i].f_send[k];
//			//ff[i_max_size][i + 1][k] = FRightRecv[i].ff_send[k];
//		}
//	}
//
//
//	//填充上下数据
//	for (int j = 0; j < MatScale[0] + 2; j++)
//	{
//		FUpSend[j].rho_send  = rho[j][j_max_size - 1];
//		FUpSend[j].u_send[0] = u  [j][j_max_size - 1][0];
//		FUpSend[j].u_send[1] = u  [j][j_max_size - 1][1];
//
//		FDownSend[j].rho_send  = rho[j][j_min];
//		FDownSend[j].u_send[0] = u  [j][j_min][0];
//		FDownSend[j].u_send[1] = u  [j][j_min][1];
//		for ( k = 0; k < Q; k++)
//		{
//			FUpSend  [j].F_send[k] = F[j][j_max_size - 1][k];
//			FUpSend[j].ff_send[k] = ff[j][j_max_size - 1][k];
//			FUpSend[j].f_send[k] = f[j][j_max_size - 1][k];
//			FDownSend[j].F_send[1] = F[j][j_min][k];
//			FDownSend[j].f_send[1] = f[j][j_min][k];
//			FDownSend[j].ff_send[1] = ff[j][j_min][k];
//		}
//	}
//
//	////数据从上往下传递，
//	//MPI_Sendrecv(FDownSend,(MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[3], 1,
//	//	         FUpRecv,  (MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[2], 1, MPI_COMM_WORLD, &Status1);
//	////数据从下往上传递
//
//	//MPI_Sendrecv(FUpSend  , (MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[2], 1,
//	//	         FDownRecv, (MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[3], 1, MPI_COMM_WORLD, &Status1);
//
//	MPI_Send(FUpSend   /*UupSend*/ ,(MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[2], 88, comm1);
//	MPI_Recv(FDownRecv/*UdownRec*/ ,(MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[3], 88, comm1, &Status1);
//					   
//	MPI_Send(FDownSend/*UdownSend*/,(MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[3], 88, comm1);
//	MPI_Recv(FUpRecv     /*UupRec*/,(MatScale[0] + 2)*sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[2], 88, comm1, &Status1);
//
//	//填充上下的数据
//	for (int j = 0; j < MatScale[0] + 2; j++)
//	{
//		rho[j][j_max_size]    = FUpRecv[j].rho_send ;
//		u  [j][j_max_size][0] = FUpRecv[j].u_send[0];
//		u  [j][j_max_size][1] = FUpRecv[j].u_send[1];
//
//		rho[j][0]    = FDownRecv[j].rho_send;
//		u  [j][0][0] = FDownRecv[j].u_send[0];
//		u  [j][0][1] = FDownRecv[j].u_send[1];
//
//		for (k = 0; k < Q; k++)
//		{
//			F[j][j_max_size][k] = FUpRecv  [j].F_send[k];
//			f[j][j_max_size][k] = FUpRecv[j].f_send[k];
//			ff[j][j_max_size][k] = FUpRecv[j].ff_send[k];
//			F[j][0][k]          = FDownRecv[j].F_send[k];
//			f[j][0][k] = FDownRecv[j].f_send[k];
//			ff[j][0][k] = FDownRecv[j].ff_send[k];
//
//		}
//
//	}
//
//}

//void    Evolution()
//{
//	for (j = j_max_size-1; j >= j_min; j--)
//	{
//		for (i = i_min; i <= i_max_size-1; i++) ////////////////////////////////////////
//		{
//			if (Soid[i][j] == 0)
//			{
//				for (int k = 0; k < Q; k++)
//					F[i][j][k] = f[i][j][k]
//					+
//					Sniu*(0.5*(f[i][j][k] + f[i][j][ne_f[k]]) - 0.5*(feq(k, rho[i][j], u[i][j]) + feq(ne_f[k], rho[i][j], u[i][j])))
//					+
//					S_q *(0.5*(f[i][j][k] - f[i][j][ne_f[k]]) - 0.5*(feq(k, rho[i][j], u[i][j]) - feq(ne_f[k], rho[i][j], u[i][j])))
//
//					+ 3.0*dt*w[k] * Force*(e_f[k][0] /** sin(theta) + e_f[k][1] * cos(theta)*/);
//			}
//		}
//	}
//	Commucate_Data();
//	MPI_Barrier(MPI_COMM_WORLD);
//	for (j = j_min; j <j_max_size; j++)     //===================advection of f 
//		for (i = i_min; i <i_max_size; i++)
//		{
//			if (Soid[i][j] == 0)
//			{
//
//				for (int k = 0; k < Q; k++)
//				{
//					ip = (i - e_f[k][0]) ;
//					jp = (j - e_f[k][1]) ;
//
//					if (Soid[ip][jp] == 2)
//					{
//
//						ff[i][j][k] = f[i][j][ne_f[k]]
//							+
//							Sniu*(0.5*(f[i][j][ne_f[k]] + f[i][j][k]) - 0.5*(feq(ne_f[k], rho[i][j], u[i][j]) + feq(k, rho[i][j], u[i][j])))
//							+
//							S_q *(0.5*(f[i][j][ne_f[k]] - f[i][j][k]) - 0.5*(feq(ne_f[k], rho[i][j], u[i][j]) - feq(k, rho[i][j], u[i][j])));
//					}
//					if (Soid[ip][jp] == 0)
//					{
//						ff[i][j][k] = F[ip][jp][k];
//					}
//				}
//			}
//
//		}
//
//
//	for (i = i_min; i <i_max_size; i++)
//		for (j = j_min; j <j_max_size; j++)
//		{
//
//			rho[i][j]  = 0.0;
//			u[i][j][0] = 0.0;
//			u[i][j][1] = 0.0;
//
//			if (Soid[i][j] == 0)
//			{
//				for (int k = 0; k < Q; k++)
//				{
//					f[i][j][k] = ff[i][j][k];
//					rho[i][j] += f[i][j][k];
//					u[i][j][0] += c*e_f[k][0] * f[i][j][k]; /////////////////////////////////////////  c!=1
//					u[i][j][1] += c*e_f[k][1] * f[i][j][k];
//					if (!(f[i][j][k] == f[i][j][k]  ))
//					{
//						cout << "the i index is :" << i << "  the j index is :" << j << " soid "<<Soid[i][j]<< endl;
//						cout << rho[i][j] << "  "<<u[i][j][0]<< endl;
//						system("pause");
//					}
//				}
//
//
//
//
//				u[i][j][0] /= rho[i][j];
//				u[i][j][1] /= rho[i][j];
//
//			}
//		}
//
//}

void    Error()   //compute error
{
	

	temp_1 = 0.0;
	temp_2 = 0.0;

	for (i = i_min; i <i_max_size; i++)
		for (j = j_min; j <j_max_size; j++)
		{
			temp_1 += sqrt((u[i][j][0] - u0[i][j][0])*(u[i][j][0] - u0[i][j][0]) + (u[i][j][1] - u0[i][j][1])*(u[i][j][1] - u0[i][j][1]));
			temp_2 += sqrt(u[i][j][0] * u[i][j][0] + u[i][j][1] * u[i][j][1]);
			u0[i][j][0] = u[i][j][0];
			u0[i][j][1] = u[i][j][1];
		}
}

void    Output_Result(int myrank)
{
	int i, j;
	ostringstream name;
	name << "size_" << Dims[0] * Dims[1] << "_rank_" << myrank << "_grid_" << NX << "_"<<coords[0]<<"_"<<coords[1]<<"_" << theta << ".dat";
	ofstream out(name.str().c_str());
	out << "Title=\"LBM Lid Driven Flow\"\n" << "VARIABLLES=\"X\",\"Y\",\"U\",\"V\"\n" << "ZONE T=\"BOX\",I=" << MatScale[0]  << ",J=" << MatScale[1]  << ",F=POINT" << endl;
	out.precision(12);
	out << setiosflags(ios::fixed);
	for (j = 1; j < j_max_size  ; j++)
	{
		for (i = 1; i < i_max_size  ; i++)
		{
			out << (double((i - 1 + Disp[0])) + 0.5) / double(NX-2) << "  " << (double((j - 1 + Disp[1])) + 0.5) / double(NY-2) << " " << u[i][j][0] << "  " << u[i][j][1];
			out << endl;
		}
		//
	}

}
void    Commucate_Data()
{
	/************************************************************************/
	/* fill the data in FLeftSend and FRightSend                            */
	/************************************************************************/
	for (int jj = j_min; jj<j_max_size; jj++)
	{
		FLeftSend[jj - 1].u_send[0] = u[i_min][jj][0];
		FLeftSend[jj - 1].u_send[1] = u[i_min][jj][1];
		FLeftSend[jj - 1].rho_send = rho[i_min][jj];
		//FLeftSend[jj - 1].soid_send = Soid[i_min][jj];



		FRightSend[jj - 1].u_send[0] = u[i_max_size - 1][jj][0];
		FRightSend[jj - 1].u_send[1] = u[i_max_size - 1][jj][1];
		FRightSend[jj - 1].rho_send = rho[i_max_size - 1][jj];
		//FRightSend[jj - 1].soid_send = Soid[i_max_size - 1][jj];


		for (int kk = 0; kk<Q; kk++)
		{
			FLeftSend[jj - 1].F_send[kk] = F[i_min][jj][kk];
			FRightSend[jj - 1].F_send[kk] = F[i_max_size - 1][jj][kk];

			//FLeftSend[jj - 1].FF_send[kk] = ff[i_min][jj][kk];
			//FRightSend[jj - 1].FF_send[kk] = ff[i_max_size - 1][jj][kk];

			//FLeftSend[jj - 1].f_send[kk] = f[i_min][jj][kk];
			//FRightSend[jj - 1].f_send[kk] = f[i_max_size - 1][jj][kk];
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);


	MPI_Sendrecv(FRightSend, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[1], 1,
		FLeftRecv, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[0], 1, comm1, &Status1);

	//数据从右向左移动
	MPI_Sendrecv(FLeftSend, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[0], 1,
		FRightRecv, MatScale[1] * sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[1], 1, comm1, &Status1);

	//左右数据填充
	for (int jj = j_min; jj<j_max_size; jj++)
	{
		u[0][jj][0] = FLeftRecv[jj - 1].u_send[0];
		u[0][jj][1] = FLeftRecv[jj - 1].u_send[1];

		u[i_max_size][jj][0] = FRightRecv[jj - 1].u_send[0];
		u[i_max_size][jj][1] = FRightRecv[jj - 1].u_send[1];

		rho[0][jj] = FLeftRecv[jj - 1].rho_send;
		rho[i_max_size][jj] = FRightRecv[jj - 1].rho_send;
		//Soid[0][jj] = FLeftRecv[jj - 1].soid_send;
		//Soid[i_max_size][jj] = FRightRecv[jj - 1].soid_send;
		for (int kk = 0; kk < Q; kk++)
		{
			//ff[0][jj][kk] = FLeftRecv[jj - 1].FF_send[kk];
			//ff[i_max_size][jj][kk] = FRightRecv[jj - 1].FF_send[kk];
			//f[0][jj][kk] = FLeftRecv[jj - 1].f_send[kk];
			//f[i_max_size][jj][kk] = FRightRecv[jj - 1].f_send[kk];

			F[0][jj][kk] = FLeftRecv[jj - 1].F_send[kk];
			F[i_max_size][jj][kk] = FRightRecv[jj - 1].F_send[kk];



		}
	}

	/************************************************************************/
	/* fill the data in FUpSend and  FDownSend                              */
	/************************************************************************/
	for (int ii = 0; ii<i_max_size + 1; ii++)
	{
		FUpSend[ii].rho_send = rho[ii][j_max_size - 1];
		//FUpSend[ii].soid_send = Soid[ii][j_max_size - 1];
		FUpSend[ii].u_send[0] = u[ii][j_max_size - 1][0];
		FUpSend[ii].u_send[1] = u[ii][j_max_size - 1][1];



		FDownSend[ii].rho_send = rho[ii][j_min];
		//FDownSend[ii].soid_send = Soid[ii][j_min];
		FDownSend[ii].u_send[0] = u[ii][j_min][0];
		FDownSend[ii].u_send[1] = u[ii][j_min][1];
		for (int k = 0; k<Q; k++)
		{
			FUpSend[ii].F_send[k] = F[ii][j_max_size - 1][k];
			FDownSend[ii].F_send[k] = F[ii][j_min][k];
			//FUpSend[ii].FF_send[k] = ff[ii][j_max_size - 1][k];
			//FDownSend[ii].FF_send[k] = ff[ii][j_min][k];

			//FUpSend[ii].f_send[k] = f[ii][j_max_size - 1][k];
			//FDownSend[ii].f_send[k] = f[ii][j_min][k];
		}
	}


	MPI_Sendrecv(FDownSend, (MatScale[0] + 2)* sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[3], 1,
		FUpRecv, (MatScale[0] + 2)* sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[2], 1, comm1, &Status1);

	//将数据从下往上移动
	MPI_Sendrecv(FUpSend, (MatScale[0] + 2)* sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[2], 1,
		FDownRecv, (MatScale[0] + 2)* sizeof(MPI_SR_Data), MPI_BYTE, Neighbor[3], 1, comm1, &Status1);

	//上下数据的填充
	MPI_Barrier(MPI_COMM_WORLD);

	for (int ii = 0; ii<i_max_size + 1; ii++)
	{
		//填充上面的数据
		u[ii][j_max_size][0] = FUpRecv[ii].u_send[0];
		u[ii][j_max_size][1] = FUpRecv[ii].u_send[1];
		//Soid[ii][j_max_size] = FUpRecv[ii].soid_send;
		rho[ii][j_max_size] = FUpRecv[ii].rho_send;
		//填充下面的数据
		u[ii][0][0] = FDownRecv[ii].u_send[0];
		u[ii][0][1] = FDownRecv[ii].u_send[1];
		//Soid[ii][0] = FDownRecv[ii].soid_send;
		rho[ii][0] = FDownRecv[ii].rho_send;

		for (int k = 0; k<Q; k++)
		{
			F[ii][j_max_size][k] = FUpRecv[ii].F_send[k];
			F[ii][0][k] = FDownRecv[ii].F_send[k];
			//ff[ii][j_max_size][k] = FUpRecv[ii].FF_send[k];
			//ff[ii][0][k] = FDownRecv[ii].FF_send[k];

			//f[ii][j_max_size][k] = FUpRecv[ii].f_send[k];
			//f[ii][0][k] = FDownRecv[ii].f_send[k];
		}
	}
}


void    Evolution()
{
	for (int i = i_min; i<i_max_size; i++)
	{
		for (int j = j_max_size - 1; j >= j_min; j--)
		{
			if (Soid[i][j] == 0)
			{
				for (int k = 0; k<Q; k++)
				{

					F[i][j][k] = f[i][j][k]
						+
						Sniu*(0.5*(f[i][j][k] + f[i][j][ne_f[k]]) - 0.5*(feq(k, rho[i][j], u[i][j]) + feq(ne_f[k], rho[i][j], u[i][j])))
						+
						S_q *(0.5*(f[i][j][k] - f[i][j][ne_f[k]]) - 0.5*(feq(k, rho[i][j], u[i][j]) - feq(ne_f[k], rho[i][j], u[i][j])))

						+ 3.0*dt*w[k] * Force*(e_f[k][0] * cos(theta) + e_f[k][1] * sin(theta));
				}
			}
		}
	}
	Commucate_Data();
	for (int i = i_min; i<i_max_size; i++)
	{
		for (int j = j_min; j<j_max_size; j++)
		{
			if (Soid[i][j] == 0)
			{
				for (int k = 0; k<Q; k++)
				{
					ip = i - e_f[k][0];
					jp = j - e_f[k][1];
					if (Soid[ip][jp] == 2)
					{
						ff[i][j][k] = f[i][j][ne_f[k]]
							+
							Sniu*(0.5*(f[i][j][ne_f[k]] + f[i][j][k]) - 0.5*(feq(ne_f[k], rho[i][j], u[i][j]) + feq(k, rho[i][j], u[i][j])))
							+
							S_q *(0.5*(f[i][j][ne_f[k]] - f[i][j][k]) - 0.5*(feq(ne_f[k], rho[i][j], u[i][j]) - feq(k, rho[i][j], u[i][j])));
					}
					if (Soid[ip][jp] == 0)
					{
						ff[i][j][k] = F[ip][jp][k];
					}
				}
			}
		}
	}

	for (int i = i_min; i<i_max_size; i++)
	{
		for (int j = j_min; j<j_max_size; j++)
		{
			rho[i][j] = 0.0;
			u[i][j][0] = 0.0;
			u[i][j][1] = 0.0;
			if (Soid[i][j] == 0)
			{
				for (int k = 0; k < Q; k++)
				{

					f[i][j][k] = ff[i][j][k];
					rho[i][j] += f[i][j][k];
					u[i][j][0] += c*e_f[k][0] * f[i][j][k]; /////////////////////////////////////////  c!=1
					u[i][j][1] += c*e_f[k][1] * f[i][j][k];
				}




				u[i][j][0] /= rho[i][j];
				u[i][j][1] /= rho[i][j];

			}

		}
	}
}
void    Drag_Force()
{

	//int myrank;
        for (j = j_min; j <j_max_size; j++)
             for (i = i_min; i <i_max_size; i++)
                     for (int k = 0; k <Q; k++)
                          F[i][j][k]=f[i][j][k];
        Commucate_Data();
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//printf("the rank %d has run j : %d  ,i : %d!\n", myrank,j_max_size,i_max_size);
	for (j = j_min; j <j_max_size; j++)     //===================advection of f 
		for (i = i_min; i <i_max_size; i++)
		{
			Drag_F[i][j][0] = 0.0;
			Drag_F[i][j][1] = 0.0;
			if (Soid[i][j] == 2)
			{
				flag = flag + 1.0;
				for (int k = 0; k < Q; k++)
				{
					ip = i + e_f[k][0];
					jp = j + e_f[k][1];
					if (Soid[ip][jp] == 0)
					{
						Drag_F[i][j][0] = Drag_F[i][j][0] + (/*F[i][j][k] +*/2.0* F[ip][jp][ne_f[k]])*double(e_f[k][0]);
						Drag_F[i][j][1] = Drag_F[i][j][1] + (/*F[i][j][k] +*/2.0* F[ip][jp][ne_f[k]])*double(e_f[k][1]);
					}
				}

			}
                        
			fx = Drag_F[i][j][0] + fx;
			fy = Drag_F[i][j][1] + fy;


		}
	/*if (fabs(flag) > 0)
	{
		printf("the rank %d has run j : %d  ,i : %d!\n", myrank, j_max_size, i_max_size);
	}*/

}

void    Soid_Flag()
{
	int *flag_left_send, *flag_right_send, *flag_up_send, *flag_down_send;
	int *flag_left_recv, *flag_right_recv, *flag_up_recv, *flag_down_recv;
	// initial the array to be sent
	flag_left_send = new int[MatScale[1]];  flag_up_send = new int[MatScale[0] + CRAT_DIM];
	flag_right_send = new int[MatScale[1]];  flag_down_send = new int[MatScale[0] + CRAT_DIM];
	flag_left_recv = new int[MatScale[1]];  flag_up_recv = new int[MatScale[0] + CRAT_DIM];
	flag_right_recv = new int[MatScale[1]];  flag_down_recv = new int[MatScale[0] + CRAT_DIM];


	int i_tid, j_tid;
	int dis_x, dis_y;
	for (int ii = i_min; ii <i_max_size; ii++)
	{
		i_tid = Disp[0] + ii - 1;
		for (int jj = j_min; jj <j_max_size; jj++)
		{
			int a = (NX - 2) / 4;
			j_tid = Disp[1] + jj - 1;
			if (i_tid<(NX - 2) / 2)
				dis_x = (NX - 2) / 2 - i_tid;
			else
				dis_x = i_tid - (NX) / 2;
			if (j_tid<(NX - 2) / 2)
				dis_y = (NY - 2) / 2 - j_tid;
			else
				dis_y = j_tid - (NY) / 2;
			if (dis_x >= (a) && dis_y >= (a))
			{
				Soid[ii][jj] = 1;
				if (dis_x == a || dis_y == a)
					Soid[ii][jj] = 2;
			}
			else
				Soid[ii][jj] = 0;
		}
	}



	//填充左右数据
	for (int i = 0; i < MatScale[1]; i++)
	{
		flag_left_send[i] = Soid[i_min][i + 1];
		flag_right_send[i] = Soid[i_max_size - 1][i + 1];
	}
	//数据从左往右传递，
	MPI_Sendrecv(flag_right_send, MatScale[1], MPI_INT, Neighbor[1], 1,
		flag_left_recv, MatScale[1], MPI_INT, Neighbor[0], 1, MPI_COMM_WORLD, &Status1);



	//数据从右往左传递，
	MPI_Sendrecv(flag_left_send, MatScale[1], MPI_INT, Neighbor[0], 1,
		flag_right_recv, MatScale[1], MPI_INT, Neighbor[1], 1, MPI_COMM_WORLD, &Status1);

	//将传递的数据进行解包
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < MatScale[1]; i++)
	{
		Soid[0][i + 1] = flag_left_recv[i];
		Soid[i_max_size][i + 1] = flag_right_recv[i];
	}


	//填充上下数据
	for (int j = 0; j < MatScale[0] + 2; j++)
	{
		flag_up_send[j] = Soid[j][j_max_size - 1];
		flag_down_send[j] = Soid[j][j_min];
	}

	//数据从上往下传递，
	MPI_Sendrecv(flag_down_send, MatScale[0] + 2, MPI_INT, Neighbor[3], 1,
		flag_up_recv, MatScale[0] + 2, MPI_INT, Neighbor[2], 1, MPI_COMM_WORLD, &Status1);
	//数据从下往上传递

	MPI_Sendrecv(flag_up_send, MatScale[0] + 2, MPI_INT, Neighbor[2], 1,
		flag_down_recv, MatScale[0] + 2, MPI_INT, Neighbor[3], 1, MPI_COMM_WORLD, &Status1);

	//填充上下的数据
	for (int j = 0; j < MatScale[0] + 2; j++)
	{
		Soid[j][j_max_size] = flag_up_recv[j];
		Soid[j][0] = flag_down_recv[j];
	}


	delete[]flag_left_send;      delete[]flag_up_send;
	delete[]flag_right_send;      delete[]flag_down_send;
	delete[]flag_left_recv;      delete[]flag_up_recv;
	delete[]flag_right_recv;      delete[]flag_down_recv;

}
