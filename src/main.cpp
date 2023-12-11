#include "grid.h"
#include "input.h"
#include <iostream>
#include <vector>
#include "chazhi.h"
#include <string>
#include <lapacke.h>
#include "scalapack_connector.h"
#include "diago.h"
#ifdef _MPI
#include "mpi.h"
#endif
#include "omp.h"
#include "timer.h"
using namespace std;
#define SLICE_SIZE(n, blk_size, my_idx, nprocs) \
    ((n) / (blk_size) / (nprocs)) +             \
        (((my_idx) < ((n) / (blk_size)) % (nprocs)) ? 1 : 0)

const char *output_filename = "output.txt";

// 输出矩阵 H 到文件中
void write_matrix_to_file(const double *H, int count)
{
    std::ofstream out_file(output_filename);
    if (out_file.is_open())
    {
        for (int i = 0; i < count; i++)
        {
            for (int j = 0; j < count; j++)
            {
                out_file << H[i * count + j] << " ";
            }
            out_file << "\n";
        }
        out_file.close();
        std::cout << "Matrix H has been written to file: " << output_filename << std::endl;
    }
    else
    {
        std::cerr << "Failed to open output file: " << output_filename << std::endl;
    }
}

int main(int argc, char *argv[])
{
    double *H_global = nullptr; // 全局矩阵数组
    string diago_lib;
    int count;
    // 初始化MPI环境
    /*#ifdef _MPI
        int ierr;
        ierr = MPI_Init(&argc, &argv);
        if (ierr != MPI_SUCCESS)
        {
            printf("MPI initialization failed!\n");
            MPI_Abort(MPI_COMM_WORLD, ierr);
        }
        // 获取进程信息
        int size, rank;
    // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);*/
    /*// Set up the grid
   int nprocs_row = std::sqrt(size);
   int nprocs_col = np / nprocs_row;
   int ctxt;
   Cblacs_pinfo(&rank, &np);
   Cblacs_get(-1, 0, &ctxt);
   Cblacs_gridinit(&ctxt, "Row-major", nprocs_row, nprocs_col);*/
    // if(rank == 0){
    // #endif
    // 根进程
    std::ofstream os;
    timer::start();
    Input input("/data/1a_bingxing/aaadazuoye/test/INPUT.txt");
    cout << "isHexahedral: " << input.isHexahedral << endl;
    cout << "lx: " << input.lx << endl;
    cout << "ly: " << input.ly << endl;
    cout << "lz: " << input.lz << endl;
    cout << "thetaxy: " << input.thetaxy << endl;
    cout << "thetayz: " << input.thetayz << endl;
    cout << "thetaxz: " << input.thetaxz << endl;
    cout << "support_SH: " << input.support_SH << endl;
    cout << "diago_lib: " << input.diago_lib << endl;
    cout << "support_Periodic_Boundary: " << input.support_Periodic_Boundary << endl;
    cout << "multi_parallel_strategies: " << input.multi_parallel_strategies << endl;
    cout << "points_path: " << input.points_path << endl;
    cout << "venergy_path: " << input.venergy_path << endl;
    cout << "distribution_path: " << input.distribution_path << endl;
    std::vector<Input::Point> points_ = input.readpointsFile(input.points_path, 50);
    count = input.count;
    diago_lib = input.diago_lib;
    // cout << "point1: (" << input.point1_[0] << ", " << input.point1_[1] << ", " << input.point1_[2] << ")" << endl;
    // cout << "point2: (" << input.point2_[0] << ", " << input.point2_[1] << ", " << input.point2_[2] << ")" << endl;

    // 读入径向分布函数
    input.readDirstibution(input.distribution_path);
    std::cout << "cutoff: " << input.cutoff << std::endl;
    double cutoff = input.cutoff;
    std::cout << "dr: " << input.dr << std::endl;
    std::cout << "mesh: " << input.mesh << std::endl;
    std::cout << "l: " << input.l << std::endl;
    std::cout << "f values: ";
    for (double val : input.f)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

     //读入V
    input.readVfile(input.venergy_path);
    std::cout << "nx : " << input.nx << std::endl;
    std::cout << "ny : " << input.ny << std::endl;
    std::cout << "nz : " << input.nz << std::endl;
    int nx = input.nx;
    int ny = input.ny;
    int nz = input.nz;
    // 避免在循环中多次访问input类的成员变量，把它存为本地变量
    double *V = new double[nx * ny * nz];
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                int index = (i * ny + j) * nz + k;
                V[index] = input.V[(i * ny + j) * nz + k];
            }
        }
    }
    cout << "成功把V存到本地啦" << endl;
    // 测试：把V都设成1
    /*
    int nx = 512;
    int ny = 512;
    int nz = 512;
    cout << "nx=" << nx << " ny = " << ny << " nz = " << nz << endl;
    double *V = new double[nx * ny * nz];
    for (int i = 0; i < nx * ny * nz; i++)
    {
        V[i] = 1;
    }*/
    int mesh = input.mesh;
    double dr = input.dr;
    // 生成距离数组d
    std::vector<double> d(mesh);
    for (int i = 0; i < mesh; i++)
    {
        d[i] = i * dr;
    }
    cout << "成功把d存到本地啦" << endl;
    double lx = input.lx;
    double ly = input.ly;
    double lz = input.lz;
    // 创建一个边长为lx ly lz，每个方向上分别有nx、ny、nz个点的网格
    Grid grid(lx, ly, lz, nx, ny, nz);
    cout << "grid初始化成功" << endl;
    // auto points = grid.generate_uniform_grid(); // 生成均匀网格
    cout << "成功生成均匀网格" << endl;
    double cell_vol = grid.cell_volume(); // 计算每个网格单元的体积
    cout << "cell_vol=" << cell_vol << endl;
    cout << "count = " << count << endl;

    Interpolator interpolator(d, input.f);
    double dx = lx / (nx - 1);
    double dy = ly / (ny - 1);
    double dz = lz / (nz - 1);

    double H[count * count] = {0}; // H矩阵

    // 尝试计算插值
    timer::tick("aaa", "jifen");
    int num_threads = 8; // 设置线程数为 8
    omp_set_num_threads(num_threads);

#pragma omp parallel for reduction(+ : H[ : count * count])
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {

                double x0 = i * dx;

                double y0 = j * dy;

                double z0 = k * dz;

                double v = V[(i * ny + j) * nz + k];//不过反正V都是1
                //double v = 1.0;
#pragma omp simd
                for (int p = 0; p < count; p++)
                {
                    // 现算p和q点的插值
                    double r1 = interpolator.distance(x0, y0, z0, points_[p].x, points_[p].y, points_[p].z);

                    if (r1 < cutoff) // 如果f1还是0就直接不算了
                    {
                        double f1 = interpolator.interpolate(r1);
#pragma omp simd
                        for (int q = 0; q <= p; q++)
                        {
                            if (q == p)
                            {
#pragma omp atomic
                                H[p * count + q] += f1 * f1 * v;
                                continue;
                            }
                            double r2 = interpolator.distance(x0, y0, z0, points_[q].x, points_[q].y, points_[q].z);

                            if (r2 < cutoff)
                            {
                                double f2 = interpolator.interpolate(r2);
#pragma omp atomic
                                H[p * count + q] += f1 * f2 * v;
                            }
                        }
                    }
                }
            }
        }
    }
    /*
    #pragma omp parallel for reduction(+ : H[ : count * count])
        for (int p = 0; p < count; p++)
        {
            for (int q = 0; q <= p; q++)
            {

                double sum = 0.0;
                if (q == p)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        double x0 = i * dx;
                        for (int j = 0; j < ny; j++)
                        {
                            double y0 = j * dy;
                            for (int k = 0; k < nz; k++)
                            {
                                double z0 = k * dz;
                                double r1 = interpolator.distance(x0, y0, z0, points_[p].x, points_[p].y, points_[p].z);
                                if (r1 < cutoff) // 如果f1还是0就直接不算了
                                {
                                    double f1 = interpolator.interpolate(r1);
                                    // double v = V[(i * ny + j) * nz + k];//不过反正V都是1
                                    double v = 1.0;
                                    sum += f1 * f1 * v;
                                }
                            }
                        }
                    }
                }
                else
                {
    #pragma omp simd
                    for (int i = 0; i < nx; i++)
                    {
                        double x0 = i * dx;
                        for (int j = 0; j < ny; j++)
                        {
                            double y0 = j * dy;
                            for (int k = 0; k < nz; k++)
                            {
                                double z0 = k * dz;
                                double r1 = interpolator.distance(x0, y0, z0, points_[p].x, points_[p].y, points_[p].z);
                                if (r1 < cutoff) // 如果f1还是0就直接不算了
                                {
                                    double f1 = interpolator.interpolate(r1);
                                    double r2 = interpolator.distance(x0, y0, z0, points_[q].x, points_[q].y, points_[q].z);
                                    if (r2 < cutoff)
                                    {
                                        double f2 = interpolator.interpolate(r2);
                                        // double v = V[(i * ny + j) * nz + k];//不过反正V都是1
                                        double v = 1.0;
                                        sum += f1 * f2 * v;
                                    }
                                }
                            }
                        }
                    }
    #pragma omp atomic
                    H[p * count + q] += sum;
                }
            }
        }*/

    delete[] V;
    // 下三角部分同上
    for (int p = 0; p < count; p++)
    {
        for (int q = 0; q < p; q++)
        {
            H[q * count + p] = H[p * count + q];
        }
    }
    for (int p = 0; p < count * count; p++)
    {
        H[p] *= cell_vol;
    }

    timer::tick("aaa", "jifen");

    // 输出初始矩阵的值
    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
            int index = i * count + j;
            std::cout << H[index] << " ";
        }
        cout << endl;
    }
    if (input.diago_lib == "lapack")
    {
        lapack_diago(H, count);
        timer::finish(os, 1);
    }

    if (input.diago_lib == "scalapack") // 分发矩阵
    {
        // Allocate and initialize the matrix H
        int n = count;
        H_global = new double[n * n];
        // 填充 H_global 数组
        for (int i = 0; i < count; i++)
        {
            for (int j = 0; j < count; j++)
            {
                H_global[i + j * count] = H[i + j * count];
            }
        }
        timer::finish(os, 1);
    }
    //}

    if (diago_lib == "scalapack")
    {
        // scalpack_diago(H_global,count,size,rank);
        //  输出矩阵 H 到文件中
        write_matrix_to_file(H_global, count);
        /*
        int n=count;
        int nb = count/np;
        int descH[9];
        int info;
            descH[0] = 1; descH[1] = 0; descH[2] = n; descH[3] = n; descH[4] = nb;
            descH[5] = nb; descH[6] = 0; descH[7] = 0; descH[8] = ctxt;
        // Diagonalize the matrix H
        double *w = new double[n];
        double *Z = new double[n*n];
        int ldz = descH[5];
        int mloc = SLICE_SIZE(count, n, rank, nprow);
        int izero = 0;    // 全局数组中第一个元素的行和列编号
        int iz = izero;                                // 矩阵 H 在全局数组中的起始列索引
        int jz = izero + ((rank - 1) * mloc);          // 矩阵 H 在全局数组中的起始行索引
        // 查询工作空间大小
        double WORK_OPTIMAL;
        pdsyev_("V", "U", &count, H_local, &IZ, &JZ, DESCH, W, Z, &IZ, &JZ, DESCH, &WORK_OPTIMAL, &LWORK, &INFO);

        // 分配工作空间
        LWORK = (int) WORK_OPTIMAL;
        double *WORK = new double[LWORK];

        pdsyev_("V", "U", &n, H_global, &izero, &izero, descH, w, Z, &izero, &izero, descH, &info);

        if (info != 0) {
            std::cerr << "Error: pdsyev_ failed with info = " << info << std::endl;
            return 1;
        }

        // Print the eigenvalues
        if (rank == 0) {
            for (int i = 0; i < n; i++) {
                std::cout << "Eigenvalue " << i << " = " << w[i] << std::endl;
            }
        }

        // Deallocate resources
        delete[] H_global;
        delete[] w;
        delete[] Z;

         // Finalize MPI
        Cblacs_gridexit(ctxt);  */
    }
    // #ifdef _MPI
    //        MPI_Finalize();
    // #endif
    return 0;
}
