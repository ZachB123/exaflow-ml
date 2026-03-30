#pragma once

#include <vector>
#include <string>

class BurgerScheme2D {
public:
    virtual ~BurgerScheme2D() = default;

    // Update both u and v for one time step
    virtual void calculateNextU(
        const std::vector<double>& u,
        const std::vector<double>& v,
        std::vector<double>& u_next,
        std::vector<double>& v_next,
        double cq,
        int Nx,
        int Ny,
        double dt,
        double dx,
        double dy,
        double kinematic_viscosity
    ) = 0;

    virtual std::string getName() const = 0;
};

class FTCS2D : public BurgerScheme2D {
public:
    void calculateNextU(
        const std::vector<double>& u,
        const std::vector<double>& v,
        std::vector<double>& u_next,
        std::vector<double>& v_next,
        double cq,
        int Nx,
        int Ny,
        double dt,
        double dx,
        double dy,
        double kinematic_viscosity
    ) override;

    std::string getName() const override;

protected:
    double calculateArtificialViscosity(
        const std::vector<double>& u,
        const std::vector<double>& v,
        double cq,
        double dx,
        double dy,
        int i,
        int j,
        int Nx,
        int Ny
    ) const;
};

class LaxWendroff2D : public BurgerScheme2D {
public:
    void calculateNextU(
        const std::vector<double>& u,
        const std::vector<double>& v,
        std::vector<double>& u_next,
        std::vector<double>& v_next,
        double cq,
        int Nx,
        int Ny,
        double dt,
        double dx,
        double dy,
        double kinematic_viscosity
    ) override;

    std::string getName() const override;

protected:
    double calculateArtificialViscosity(
        const std::vector<double>& u,
        const std::vector<double>& v,
        double cq,
        double dx,
        double dy,
        int i,
        int j,
        int Nx,
        int Ny
    ) const;
};

class Godunov2D : public BurgerScheme2D {
public:
    void calculateNextU(
        const std::vector<double>& u,
        const std::vector<double>& v,
        std::vector<double>& u_next,
        std::vector<double>& v_next,
        double cq,
        int Nx,
        int Ny,
        double dt,
        double dx,
        double dy,
        double kinematic_viscosity
    ) override;

    std::string getName() const override;

protected:
    double godunovFlux(double qLeft, double qRight) const;
    double transportFlux(double wLeft, double wRight,
                         double qLeft, double qRight) const;
};