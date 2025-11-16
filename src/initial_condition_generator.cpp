#include <cmath>
#include <sstream>
#include <iomanip>

#include "initial_condition_generator.h"


RandomInitialCondition::RandomInitialCondition(int n, double domain_length, bool always_positive, bool wrapAround, unsigned seed)
    : domain_length_(domain_length), bias_(0.0)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> amplitude_distribution(RandomInitialCondition::AMP_MIN, RandomInitialCondition::AMP_MAX);
    std::uniform_real_distribution<double> phase_shift_distribution(0.0, domain_length);

    terms_.reserve(n);

    double vertical_shift = 0.0;

    for (int i = 0; i < n; i++) {

        double current_amplitude = amplitude_distribution(rng);
        double current_phase_shift = phase_shift_distribution(rng);
        double current_frequency;

        if (wrapAround) {
            // Fourier frequency: b = 2πk / L
            std::uniform_int_distribution<int> integral_frequency_distribution(RandomInitialCondition::WRAP_AROUND_FREQUENCY_MULTIPLIER_MIN, RandomInitialCondition::WRAP_AROUND_FREQUENCY_MULTIPLIER_MAX);
            int k = integral_frequency_distribution(rng);
            current_frequency = 2.0 * M_PI * k / domain_length_;
        } else {
            // random frequency in a normal range
            std::uniform_real_distribution<double> frequency_distribution(RandomInitialCondition::FREQUENCY_MULTIPLIER_MIN, RandomInitialCondition::FREQUENCY_MULTIPLIER_MAX);
            current_frequency = frequency_distribution(rng);
        }

        terms_.push_back({current_amplitude, current_frequency, current_phase_shift});
        vertical_shift += std::abs(current_amplitude);
    }

    if (always_positive) {
        // slight overshoot factor to ensure positivity
        // only set the bias if always positive, otherwise we will just add 0
        bias_ = vertical_shift * 1.1;
    }
}

double RandomInitialCondition::operator()(double x) const {
    double sum = bias_;
    for (const auto& t : terms_) {
        sum += t.amplitude * std::sin(t.frequency * (x - t.phase_shift));
    }
    return sum;
}

std::string RandomInitialCondition::toString() const {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);

    out << "f(x) = ";

    // Bias term (only if nonzero)
    if (std::abs(bias_) > 1e-12) {
        out << bias_;
    } else if (!terms_.empty()) {
        out << "0";
    }

    // Add sinusoid terms
    for (const auto& t : terms_) {
        out << " + "
            << t.amplitude << "*sin("
            << t.frequency << "*(x - "
            << t.phase_shift << "))";
    }

    return out.str();
}