#ifndef SCALAPACK_DIAGO_H
#define SCALAPACK_DIAGO_H
#include "scalapack_connector.h"
#include "timer.h"
// 引入 MPI 头文件
#include <mpi.h>
void lapack_diago(double*H,int count){
    timer::tick("aaaa", "lapack");
    // 调用LAPACK库函数dsyev计算特征值和特征向量
        char jobz = 'V'; // 计算特征向量
        char uplo = 'U'; // 上三角部分存储矩阵H
        int n = count;   // 矩阵H维度为2
        double w[count]; // 存储特征值
        int lda = n;

        // 执行对角化
        int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, n, H, lda, w);

        if (info != 0)
        {
            cout << "Failed to diagonalize matrix." << endl;
        }
        // 输出特征值和特征向量
        for (int i = 0; i < n; i++) {
            std::cout << "Eigenvalue " << i << " = " << w[i] << std::endl;
        }
        //std::cout << "特征值为: " << w[0] << ", " << w[1] << std::endl;
        // std::cout << "第一个特征向量为: [" << vl[0] << ", " << vl[1] << "]" << std::endl;
        // std::cout << "第二个特征向量为: [" << vl[2] << ", " << vl[3] << "]" << std::endl;
        std::cout << "特征向量为：" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int index = i * n + j;
                std::cout << H[index] << " ";
            }
            cout << endl;
        }
        timer::tick("aaaa", "lapack");
}
void scalpack_diago(double *H,int m,int mpi_size,int mpi_rank){
    MPI_Status status;
    MPI_Bcast(H,m*m, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    int ictxt,nprow=2,npcol=2,myrow,mycol,nb=1;//
    int info,itemp;
    int izero=0, ione=1;
	char Row[4]="Row",V='V',U='U';
    
    Cblacs_pinfo(&mpi_rank, &mpi_size);//tell scalapack size and rank
    Cblacs_get(-1, 0, &ictxt);//initialize ictxt
    Cblacs_gridinit(&ictxt, Row, nprow, npcol);//tell scalapack nprow and npcol
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);//get my row and my col
    int descA[9],descZ[9];
    int Anrow = numroc_(&m,&nb,&myrow,&izero,&nprow);//get row width of this grid
    int Ancol = numroc_(&m,&nb,&mycol,&izero,&npcol);//get col width of this grid
    int Znrow=Anrow,Zncol=Ancol;

    descinit_(descA,&m,&m,&nb,&nb,&izero,&izero,&ictxt,&Anrow,&info);//initialize descA which store the messages of A
    descinit_(descZ,&m,&m,&nb,&nb,&izero,&izero,&ictxt,&Znrow,&info);//initialize descZ which store the messages of Z
    double *A = new double[Anrow*Ancol];
    double *Z = new double[Znrow*Zncol];
    
    int* Anrow_array = new int[mpi_size];
    int* Ancol_array = new int[mpi_size];
    int* myrow_array = new int[mpi_size];
    int* mycol_array = new int[mpi_size];

    MPI_Allgather(&Anrow, 1, MPI_INT, Anrow_array, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&Ancol, 1, MPI_INT, Ancol_array, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&myrow, 1, MPI_INT, myrow_array, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&mycol, 1, MPI_INT, mycol_array, 1, MPI_INT, MPI_COMM_WORLD);

    int *flag1=new int[nprow];
	int *flag2=new int[npcol];
    int start_position[2]={};
    for(int i=0;i<nprow;i++)flag1[i]=0;
    for(int i=0;i<npcol;i++)flag2[i]=0;
    for(int i=0;i<mpi_size;i++)
    	if(myrow_array[i]<myrow&&flag1[myrow_array[i]]==0)
		{
        	start_position[0]+=Anrow_array[i];
        	flag1[myrow_array[i]]=1;
    	}
	for(int i=0;i<mpi_size;i++)
    	if(mycol_array[i]<mycol&&flag2[mycol_array[i]]==0)
		{
        	start_position[1]+=Ancol_array[i];
        	flag2[mycol_array[i]]=1;
    	}
    for(int i=0;i<Anrow;i++) 
        for(int j=0;j<Ancol;j++) 
		{
            int y_in_H=start_position[0]+i;
            int x_in_H=start_position[1]+j;
            A[i*Ancol+j]=H[y_in_H*m+x_in_H];
        }

    double *work=new double[1];
    int lwork = -1;//set lwork = -1 to query work space
    double *w = new double[m];
    //double *Z = new double[m*m];
    pdsyev_(&V,&U,&m,A,&ione,&ione,descA,w,Z,&ione,&ione,descZ,work,&lwork,&info);
    lwork = (int)work[0];
    delete[] work;
    work = new double[lwork];
    pdsyev_(&V,&U,&m,A,&ione,&ione,descA,w,Z,&ione,&ione,descZ,work,&lwork,&info);//call pdsyev to diagonize it 
    double *TEMP = new double[m*m];
    for(int i=0;i<m*m;i++)
	{
		H[i]=0.0;
		TEMP[i]=0.0;
	}
	for(int i=0;i<Znrow;i++) 
        for(int j=0;j<Zncol;j++) 
		{
            int y_in_H=start_position[0]+i;
            int x_in_H=start_position[1]+j;
            H[y_in_H*m+x_in_H]=Z[i*Zncol+j];
        }
    MPI_Barrier(MPI_COMM_WORLD);
	for(int i=1;i<mpi_size;i++)
	{
		if(i==mpi_rank)
		{
			MPI_Send(H,m*m,MPI_DOUBLE,0,92,MPI_COMM_WORLD);
		}
		if(mpi_rank==0)
		{
			MPI_Recv(TEMP,m*m,MPI_DOUBLE,i,92,MPI_COMM_WORLD,&status);
			for(int i=0;i<m*m;i++)
			{
				H[i]+=TEMP[i];
			}
			
		}
	}
    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi_rank == 0){
        for (int i = 0; i < m; i++) {
            std::cout << "Eigenvalue " << i << " = " << w[i] << std::endl;
        }
    }
	delete[] A;
	delete[] Z;
    delete[] TEMP;
	delete[] Anrow_array;
	delete[] Ancol_array;
	delete[] myrow_array;
	delete[] mycol_array;
	delete[] flag1;
	delete[] flag2;
	delete[] work;
	A=nullptr;
	Z=nullptr;
	Anrow_array=nullptr;
	Ancol_array=nullptr;
	myrow_array=nullptr;
	mycol_array=nullptr;
	flag1=nullptr;
	flag2=nullptr;
	work=nullptr;

}
#endif