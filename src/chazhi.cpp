#include "chazhi.h"
#include "timer.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

double Interpolator::distance(double x, double y, double z, double cx, double cy, double cz) {
        double dx = x - cx;
        double dy = y - cy;
        double dz = z - cz;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

 double Interpolator::interpolate(double r) {
        timer::tick("Interpolator","interpolate");
        int n = r_.size();
        int i = 1;
        while (i < n && r_[i] < r) i++;
        double h = r_[i] - r_[i-1];
        double t = (r - r_[i-1]) / h;
        double t2 = t * t;
        double t3 = t2 * t;
        double a = 2.0 * t3 - 3.0 * t2 + 1.0;
        double b = t3 - 2.0 * t2 + t;
        double c = -2.0 * t3 + 3.0 * t2;
        double d = t3 - t2;
        return a * f_[i-1] + b * h * f_[i] + c * f_[i+1] + d * h * f_[i+2];
        timer::tick("Interpolator","interpolate");
    }
