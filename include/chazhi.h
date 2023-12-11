#ifndef CHAZHI_H
#define CHAZHI_H

#include <vector>
#include <math.h>

#include <vector>
#include <cmath>

class Interpolator {
public:
    std::vector<double> r_;
    std::vector<double> f_;

    Interpolator(const std::vector<double>& r, const std::vector<double>& f) : 
        r_(r), f_(f) {}

    double distance(double x, double y, double z, double cx, double cy, double cz) ;

    double interpolate(double r);
    

};
#endif