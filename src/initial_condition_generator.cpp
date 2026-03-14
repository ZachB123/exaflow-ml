#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>


#include "initial_condition_generator.h"
#include "math_utils.h"


RandomInitialCondition::RandomInitialCondition(const RandomInitialConditionConfig& config, bool always_positive, bool wrapAround, unsigned seed)
    : 
        config(config),
        domain_length(config.domain_length), 
        bias(0.0)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> amplitude_distribution(config.amp_min, config.amp_max);
    std::uniform_real_distribution<double> phase_shift_distribution(0.0, domain_length);

    terms.reserve(config.n);

    double vertical_shift = 0.0;

    for (int i = 0; i < config.n; i++) {
        double current_amplitude = amplitude_distribution(rng);
        double current_phase_shift = phase_shift_distribution(rng);
        double current_frequency;

        if (wrapAround) {
            // Fourier frequency: b = 2πk / L
            std::uniform_int_distribution<int> integral_frequency_distribution(config.wrap_around_frequency_multiplier_min, config.wrap_around_frequency_multiplier_max);
            int k = integral_frequency_distribution(rng);
            current_frequency = 2.0 * M_PI * k / domain_length;
        } else {
            // random frequency in a normal range
            std::uniform_real_distribution<double> frequency_distribution(config.frequency_multiplier_min, config.frequency_multiplier_max);
            current_frequency = frequency_distribution(rng);
        }

        terms.push_back({current_amplitude, current_frequency, current_phase_shift});
        vertical_shift += std::abs(current_amplitude);
    }

    if (always_positive) {
        // slight overshoot factor to ensure positivity
        // only set the bias if always positive, otherwise we will just add 0
        bias = vertical_shift * 1.1;
    }

    max_u = compute_max_abs(*this, domain_length);
}

double RandomInitialCondition::operator()(double x) const {
    double sum = bias;
    for (const auto& t : terms) {
        sum += t.amplitude * std::sin(t.frequency * (x - t.phase_shift));
    }
    return sum;
}

double RandomInitialCondition::getMaxU() const {
    return max_u;
}

void RandomInitialCondition::appendMetadata(nlohmann::json& metadata) const {
    metadata["config"] = {
        {"n", config.n},
        {"domain_length", config.domain_length},
        {"amp_min", config.amp_min},
        {"amp_max", config.amp_max},
        {"frequency_multiplier_min", config.frequency_multiplier_min},
        {"frequency_multiplier_max", config.frequency_multiplier_max},
        {"wrap_around_frequency_multiplier_min", config.wrap_around_frequency_multiplier_min},
        {"wrap_around_frequency_multiplier_max", config.wrap_around_frequency_multiplier_max}
    };
    nlohmann::json terms_array = nlohmann::json::array();
    for (const Term& term : terms) {
        terms_array.push_back({
            {"amplitude", term.amplitude},
            {"frequency", term.frequency},
            {"phase_shift", term.phase_shift}
        });
    }
    metadata["terms"] = terms_array;
    metadata["bias"] = bias;
    metadata["max_u"] = max_u;
}

std::string RandomInitialCondition::toString() const {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);

    out << "f(x) = ";

    // Bias term (only if nonzero)
    if (std::abs(bias) > 1e-12) {
        out << bias;
    } else if (!terms.empty()) {
        out << "0";
    }

    // Add sinusoid terms
    for (const auto& t : terms) {
        out << " + "
            << t.amplitude << "*sin("
            << t.frequency << "*(x - "
            << t.phase_shift << "))";
    }

    return out.str();
}