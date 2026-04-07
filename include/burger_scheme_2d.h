#pragma once

#include <vector>
#include <string>

class BurgerScheme2D {
    public:
        virtual ~BurgerScheme2D() = default;

        // Update both u and v for one time step
        virtual void calculateNextUandV(
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

    protected:
        // function to calculate artificial viscosity, implemented in subclasses
        virtual double calculateArtificialViscosity(
            const std::vector<double>&, // u
            const std::vector<double>&, // v
            double, // cq 
            double, // spatial_step_size_x
            double, // spatial_step_size_y
            int, // i
            int, // j
            int, // num_domain_points_x
            int ) const // num_domain_points_y
        {
            return 0.0; // default: no artificial viscosity
        }
};

class FTCS2D : public BurgerScheme2D {
    public:
        void calculateNextUandV(
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
        ) const override;
};

class LaxWendroff2D : public BurgerScheme2D {
    public:
        void calculateNextUandV(
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
            double spatial_step_size_x,
            double spatial_step_size_y,
            int i,
            int j,
            int num_domain_points_x,
            int num_domain_points_y) const override;
};

class Godunov2D : public BurgerScheme2D {
    public:
        void calculateNextUandV(
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