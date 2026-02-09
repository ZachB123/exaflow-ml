#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <limits>


#include "initial_condition_generator.h"


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
}

double RandomInitialCondition::operator()(double x) const {
    double sum = bias;
    for (const auto& t : terms) {
        sum += t.amplitude * std::sin(t.frequency * (x - t.phase_shift));
    }
    return sum;
}


void RandomInitialCondition::saveMetadataJSON(const std::filesystem::path& base_path,
                                              const std::filesystem::path& sample_folder,
                                              int time_steps,
                                              double time_step_size,
                                              int num_domain_points,
                                              double spatial_step_size,
                                              const std::string& stencil_name) const
{
    std::filesystem::path folder = base_path / sample_folder;

    // Create directory tree if not exists
    std::filesystem::create_directories(folder);

    std::filesystem::path filepath = folder / "metadata.json";

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open metadata file: " << filepath << std::endl;
        return;
    }

    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    file << "{\n";
    file << "  \"config\": {\n";
    file << "    \"n\": " << config.n << ",\n";
    file << "    \"domain_length\": " << config.domain_length << ",\n";
    file << "    \"amp_min\": " << config.amp_min << ",\n";
    file << "    \"amp_max\": " << config.amp_max << ",\n";
    file << "    \"frequency_multiplier_min\": " << config.frequency_multiplier_min << ",\n";
    file << "    \"frequency_multiplier_max\": " << config.frequency_multiplier_max << ",\n";
    file << "    \"wrap_around_frequency_multiplier_min\": " << config.wrap_around_frequency_multiplier_min << ",\n";
    file << "    \"wrap_around_frequency_multiplier_max\": " << config.wrap_around_frequency_multiplier_max << "\n";
    file << "  },\n\n";

    file << "  \"terms\": [\n";
    for (size_t i = 0; i < terms.size(); ++i) {
        const auto& t = terms[i];
        file << "    {\n";
        file << "      \"amplitude\": " << t.amplitude << ",\n";
        file << "      \"frequency\": " << t.frequency << ",\n";
        file << "      \"phase_shift\": " << t.phase_shift << "\n";
        if (i + 1 < terms.size())
            file << "    },\n";
        else
            file << "    }\n";
    }
    file << "  ],\n\n";

    file << "  \"bias\": " << bias << ",\n\n";

    file << "  \"solver\": {\n";
    file << "    \"time_steps\": " << time_steps << ",\n";
    file << "    \"time_step_size\": " << time_step_size << ",\n";
    file << "    \"num_domain_points\": " << num_domain_points << ",\n";
    file << "    \"spatial_step_size\": " << spatial_step_size << ",\n";
    file << "    \"stencil_name\": \"" << stencil_name << "\"\n";
    file << "  }\n";
    file << "}\n";
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