#pragma once

#include <vector>
#include <string>

class BurgerStencil {
public:

    // destructor
    virtual ~BurgerStencil() = default;

    // = 0 enforces its implementation in subclasses
    virtual void calculateNextU(
                const std::vector<double>& u, 
                std::vector<double>& u_next, 
                double cq, 
                int num_domain_points, 
                double time_step_size, 
                double spatial_step_size, 
                double kinematic_viscosity) = 0;

    // get the name of the stencil
    virtual std::string getName() const = 0;

protected:
    virtual double calculateArtificialViscosity(
        const std::vector<double>&, // u
        double, // cq
        double, // spatial_step_size
        int, //i
        int ) const //num_domain_points
    {
        return 0.0;   // default: no artificial viscosity
    }
};

class FTCS : public BurgerStencil {
public:
    void calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) override;

    std::string getName() const override;

protected:
    double calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i, int num_domain_points) const override;
};

class LaxWendroff : public BurgerStencil {
public:
    void calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) override;

    std::string getName() const override;

protected:
    double calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i, int num_domain_points) const override;
};

class Godunov : public BurgerStencil {
public:
    void calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double /*cq (unused)*/, int n, double dt, double dx, double kinematic_viscosity) override;
    
    std::string getName() const override;
    
protected:
    double godunovFlux(double uL, double uR) const;
};
