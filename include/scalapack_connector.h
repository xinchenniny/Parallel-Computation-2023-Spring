#ifndef SCALAPACK_CONNECTOR_H
#define SCALAPACK_CONNECTOR_H
// 引入 MPI 头文件
#include <mpi.h>
extern "C"
{
    void pdsyev_(const char *jobz, const char *uplo, int *n, double *a,
                 int *ia, int *ja, int *desca, double *w, double *z, int *iz,
                 int *jz, int *descz, double *work, int *lwork, int *info);
    void descinit_(int *desc, int *m, int *n, int *mb, int *nb,
                   int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);
    void Cblacs_gridinit(int *Context, const char *layout, int nprow, int npcol);
    void Cblacs_pcoord(int context, int proc, int *iproc, int *jproc);
    void Cblacs_gridexit(int context);
    void Cblacs_get(int context, int request, int *value);
    void Cblacs_pinfo(int *mypnum, int *nprocs);
    void Cblacs_gridinfo(int context, int *nprow, int *npcol, int *myrow, int *mycol);
    int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
    void pdgemr2d_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
                   double *b, int *ib, int *jb, int *descb,
                   int *context, int *info);
    int Cblacs_exit(int)
    {
        int ictxt;
        Cblacs_get(-1, 0, &ictxt);
        MPI_Finalize();
        return 0;
    }
}
#endif