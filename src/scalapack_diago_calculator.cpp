#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include "../include/scalapack_connector.h"

const char *input_filename = "../output.txt";
int count = 0;

// 从文件中读取矩阵 H 的值
void read_matrix_from_file(double *H_global)
{
    std::ifstream in_file(input_filename);
    if (in_file.is_open())
    {
        double value;
        int i = 0, j = 0;
        while (in_file >> value)
        {
            H_global[i * count + j] = value;
            if (j == count)
            {
                j = 0;
                i++;
            }else{
                j++;

            }
        }
        in_file.close();
        std::cout << "Matrix H has been read from file: " << input_filename << std::endl;
    }
    else
    {
        std::cerr << "Failed to open input file: " << input_filename << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv)
{
    // 初始化 MPI
    MPI_Init(&argc, &argv);

    // 获取当前进程的排名和总数
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // 定义变量
    int info, nprow = num_procs, npcol = num_procs, nb = 25, lda, lwork, liwork;
    nprow = 2; // 进程网格中的行数
    npcol = 2; // 进程网格中的列数

    double d_one = 1.0, d_zero = 0.0;
    double *A_global = nullptr, *A_local = nullptr, *W = nullptr, *Z = nullptr, *work = nullptr;
    int *iwork = nullptr, descA[9], descZ[9], ctxt;
    int izero = 1;

    // 确定矩阵的维度和数组长度
    count = 50;
    lda = count / nprow;
    if (count % nprow != 0)
        lda++;

    // 在当前进程上分配 A_global 数组
    if (rank == 0)
    {
        A_global = new double[count * count];
        read_matrix_from_file(A_global);
    }

    // 分发数据到各个进程
    A_local = new double[lda * count+1000];
    std::fill_n(A_local, lda * count, 0.0);
    Cblacs_pinfo(&rank, &num_procs);
    Cblacs_get(0, 0, &ctxt);
    Cblacs_gridinit(&ctxt, "Row-major", nprow, npcol);
    int i,j;
    Cblacs_gridinfo(ctxt, &nprow, &npcol, &i, &j);//get my row and my col

    descinit_(descA, &count, &count, &nb, &nb, &izero, &izero, &ctxt, &lda, &info);

    int descB[9];
    int myrow, mycol;
    myrow = numroc_(&count, &nb, &i, &izero, &nprow);
    mycol = numroc_(&count, &nb, &j, &izero, &npcol);
    std::cout << myrow <<" " <<mycol<<std::endl; 

    descinit_(descB, &count, &count, &nb, &nb, &izero, &izero, &ctxt, &myrow, &info);

    //这个函数搞死我了……
        pdgemr2d_(&nb, &nb, A_global, &izero, &izero, descA, A_local, &myrow, &mycol, descB, &ctxt, &info);
    // 定义工作区和描述符
    W = new double[count];
    Z = new double[lda * count];
    lwork = 4 * count;
    work = new double[lwork];
    liwork = 4 * count;
    iwork = new int[liwork];
    descinit_(descZ, &count, &count, &nb, &nb, &izero, &izero, &ctxt, &lda, &info);

    // 进行矩阵对角化
    const char jobz = 'V';
    const char uplo = 'U';
    pdsyev_(&jobz, &uplo, &count, A_local, &izero, &izero, descA, W, Z, &izero, &izero, descZ, work, &lwork, &info);

    // 输出特征值和特征向量
    if (rank == 0)
    {
        std::cout << "Eigenvalues: ";
        for (int i = 0; i < count; i++)
        {
            std::cout << W[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Eigenvectors: " << std::endl;
        for (int i = 0; i < count; i++)
        {
            for (int j = 0; j < count; j++)
            {
                std::cout << Z[i * lda + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // 释放资源并退出 MPI
    delete[] A_global;
    delete[] A_local;
    delete[] W;
    delete[] Z;
    delete[] work;
    delete[] iwork;
    Cblacs_gridexit(ctxt);
    Cblacs_exit(0);
    MPI_Finalize();
    return 0;
}