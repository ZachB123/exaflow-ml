#include "math_utils.h"
#include <cmath>
#include <algorithm>

double compute_max_abs(const std::function<double(double)>& f, double domain_length, int num_samples) {
    double max_val = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        double x = domain_length * i / (num_samples - 1);
        max_val = std::max(max_val, std::abs(f(x)));
    }
    return max_val;
}
