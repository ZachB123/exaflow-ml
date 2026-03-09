#include <iostream>
#include <vector>
#include "burger_scheme_2d.h"

//'regular function instead for now, can redefine when FTCS2D is a class'
void calculateNextU(const std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& u_next, double cq, int nx, int ny, double dt, double dx, double dy, double nu){
    