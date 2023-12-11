#include "grid.h"
#include <cmath>
#include <iostream>
#include "timer.h"
using namespace std;

Grid::Grid(double lx, double ly, double lz, int nx, int ny, int nz)
    : m_lx(lx), m_ly(ly), m_lz(lz), m_nx(nx), m_ny(ny), m_nz(nz) {}

vector<vector<double>> Grid::generate_uniform_grid() {
    timer::tick("Grid","generate_uniform_grid");
    vector<vector<double>> points;
    double dx = m_lx / (m_nx-1);
    double dy = m_ly / (m_ny-1);
    double dz = m_lz / (m_nz-1);

    for (int iz = 0; iz < m_nz; ++iz) {
        for (int iy = 0; iy < m_ny; ++iy) {
            for (int ix = 0; ix < m_nx; ++ix) {
                double x = ix * dx;
                double y = iy * dy;
                double z = iz * dz;
                points.push_back({x, y, z});
            }
        }
    }
    timer::tick("Grid","generate_uniform_grid");
    return points;
}

double Grid::cell_volume() {
    timer::tick("Grid","cell_volume");
    double dx = m_lx / (m_nx-1);
    double dy = m_ly / (m_ny-1);
    double dz = m_lz / (m_nz-1);
    timer::tick("Grid","cell_volume");
    return dx * dy * dz;
}
