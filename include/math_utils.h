#pragma once
#include <functional>

double compute_max_abs(const std::function<double(double)>& f, double domain_length, int num_samples = 10000);
