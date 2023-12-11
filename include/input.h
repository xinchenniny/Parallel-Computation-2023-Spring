#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>


using namespace std;

class Input {
public:
    Input(string filename);

    bool isHexahedral;
    double lx;
    double ly;
    double lz;
    double thetaxy;
    double thetayz;
    double thetaxz;
    bool support_SH;
    string diago_lib;
    int support_Periodic_Boundary;
    int multi_parallel_strategies;
    string points_path;
    string venergy_path;
    string distribution_path;
    vector<double> point1_;
    vector<double> point2_;
    struct Point {
    double x;
    double y;
    double z;
    };
    int count;//点坐标的个数
    std::vector<Point> readpointsFile(const std::string& filename, int limit);
    //void readpointsFile(const std::string& filename);
    double cutoff;
    double dr;
    int mesh;
    int l;
    std::vector<double> f;
    int nx,ny,nz;
    void readDirstibution(const std::string& filename);
    std::vector<double> V;

    void readVfile(const std::string& filename);
private:
    void parseFile(string filename);
};

#endif // INPUT_H
