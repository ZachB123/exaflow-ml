#pragma once

#include <vector>
#include <string>

class BurgerScheme2D {
    public:
        // destructor
        virtual ~BurgerScheme2D() = default;

        // function for calculateNextU, implemented in subclasses
        // = 0 enforces its implementation in the subclasses
        virtual void calculateNextU(
                    const std::vector<std::vector<double>>& u, 
                    std::vector<std::vector<double>>& u_next, 
                    double cq, 
                    int num_domain_points_x, 
                    int num_domain_points_y, 
                    double time_step_size,
                    double spatial_step_size_x, 
                    double spatial_step_size_y, 
                    double kinematic_viscosity) = 0;

        // get the name of the scheme
        virtual std::string getName() const = 0;

    protected:
        // function to calculate artifiical viscosity, implemented in subclasses
        virtual double calculateArtificialViscosity(
            const std::vector<std::vector<double>>&, // u 
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
        void calculateNextU(
                    const std::vector<std::vector<double>>& u, 
                    std::vector<std::vector<double>>& u_next, 
                    double cq, 
                    int num_domain_points_x, 
                    int num_domain_points_y, 
                    double time_step_size,
                    double spatial_step_size_x, 
                    double spatial_step_size_y, 
                    double kinematic_viscosity) override;

        std::string getName() const override;
    
    protected:
        double calculateArtificialViscosity(
            const std::vector<std::vector<double>>& u, 
            double cq,
            double spatial_step_size_x,
            double spatial_step_size_y,
            int i,
            int j,
            int num_domain_points_x,
            int num_domain_points_y) const override;
};

class LaxWendroff2D : public BurgerScheme2D {
    public:
        void calculateNextU(
                    const std::vector<std::vector<double>>& u, 
                    std::vector<std::vector<double>>& u_next, 
                    double cq, 
                    int num_domain_points_x, 
                    int num_domain_points_y, 
                    double time_step_size,
                    double spatial_step_size_x, 
                    double spatial_step_size_y, 
                    double kinematic_viscosity) override;

        std::string getName() const override;
    
    protected:
        double calculateArtificialViscosity(
            const std::vector<std::vector<double>>& u,
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
        void calculateNextU(
                    const std::vector<std::vector<double>>& u, 
                    std::vector<std::vector<double>>& u_next, 
                    // double cq, (unused)
                    int num_domain_points_x, 
                    int num_domain_points_y, 
                    double time_step_size,
                    double spatial_step_size_x, 
                    double spatial_step_size_y, 
                    double kinematic_viscosity) override;

        std::string getName() const override;

    protected:
        double godunovFlux(double u_left, double r_right, double u_up, double u_down) const;
};