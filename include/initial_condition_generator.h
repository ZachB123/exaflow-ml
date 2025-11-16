#pragma once
#include <vector>
#include <random>
#include <ostream>

class RandomInitialCondition {

public:
    static constexpr double AMP_MIN = -1.0;
    static constexpr double AMP_MAX =  1.0;

    static constexpr int WRAP_AROUND_FREQUENCY_MULTIPLIER_MIN = 1;
    static constexpr int WRAP_AROUND_FREQUENCY_MULTIPLIER_MAX = 5;

    static constexpr double FREQUENCY_MULTIPLIER_MIN = 0.1;
    static constexpr double FREQUENCY_MULTIPLIER_MAX = 5.0;

    struct Term {
        double amplitude;
        double frequency;
        double phase_shift;
    };

    RandomInitialCondition(int n,
                   double domainLength,
                   bool alwaysPositive = false,
                   bool wrapAround = false,
                   unsigned seed = std::random_device{}());

    double operator()(double x) const;

    std::string toString() const;

private:
    double domain_length_;
    std::vector<Term> terms_;
    double bias_; // used only if alwaysPositive = true
};
