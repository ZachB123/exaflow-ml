#pragma once

#include <vector>

class BurgerStencil {
public:

    // destructor
    virtual ~BurgerStencil() = default;

    // = 0 enforces its implementation in subclasses
    virtual void calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) = 0;

protected:
    virtual double calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i) const = 0;
};

class FTCS : public BurgerStencil {
public:
    void calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) override;

protected:
    double calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i) const override;
};

class LaxWendroff : public BurgerStencil {
public:
    void calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) override;

protected:
    double calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i) const override;
};