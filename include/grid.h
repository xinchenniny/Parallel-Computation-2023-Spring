#ifndef GRID_H
#define GRID_H

#include <vector>

class Grid {
public:
    // 构造函数，输入网格边长 lx, ly, lz 和每个方向点数 nx, ny, nz
    Grid(double lx, double ly, double lz, int nx, int ny, int nz);

    // 生成均匀网格
    std::vector<std::vector<double>> generate_uniform_grid();

    // 计算每个网格单元的体积
    double cell_volume();

private:
    double m_lx; // 网格边长 lx
    double m_ly; // 网格边长 ly
    double m_lz; // 网格边长 lz
    int m_nx;    // 每个方向点数 nx
    int m_ny;    // 每个方向点数 ny
    int m_nz;    // 每个方向点数 nz
};

#endif // GRID_H
