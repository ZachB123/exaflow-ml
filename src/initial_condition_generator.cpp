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
    int /*time_steps*/,
    double /*time_step_size*/,
    int /*num_domain_points*/,
    double /*spatial_step_size*/,
    const std::string& /*scheme_name*/) const
{
    std::filesystem::path filepath = (base_path / sample_folder) / "metadata.json";

    // Read existing solver metadata
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "Failed to open solver metadata for reading: " << filepath << "\n";
        return;
    }
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    // Find last non-whitespace char and ensure it's a closing brace
    auto is_whitespace = [](unsigned char c) { return std::isspace(c) != 0; };// helper function to check for whitespace
    size_t end = content.size();
    while (end > 0 && is_whitespace(static_cast<unsigned char>(content[end - 1]))) --end;
    if (end == 0 || content[end - 1] != '}') {
        std::cerr << "Invalid solver metadata JSON (missing closing '}'): " << filepath << "\n";
        return;
    }

    // Skip if data already appended
    if (content.find("\"config\"") != std::string::npos ||
        content.find("\"terms\"") != std::string::npos ||
        content.find("\"bias\"") != std::string::npos) {
        std::cerr << "IC metadata already present. Skipping: " << filepath << "\n";
        return;
    }

    // Remove final '}' and trailing whitespace before it
    content.erase(end - 1);
    while (!content.empty() && is_whitespace(static_cast<unsigned char>(content.back()))) content.pop_back();

    // Build the append fragment
    std::ostringstream a;
    a << std::setprecision(std::numeric_limits<double>::max_digits10);

    a << ",\n\n" << "  \"config\": {\n";
    a << "    \"n\": " << config.n << ",\n";
    a << "    \"domain_length\": " << config.domain_length << ",\n";
    a << "    \"amp_min\": " << config.amp_min << ",\n";
    a << "    \"amp_max\": " << config.amp_max << ",\n";
    a << "    \"frequency_multiplier_min\": " << config.frequency_multiplier_min << ",\n";
    a << "    \"frequency_multiplier_max\": " << config.frequency_multiplier_max << ",\n";
    a << "    \"wrap_around_frequency_multiplier_min\": " << config.wrap_around_frequency_multiplier_min << ",\n";
    a << "    \"wrap_around_frequency_multiplier_max\": " << config.wrap_around_frequency_multiplier_max << "\n";
    a << "  },\n\n";

    a << "  \"terms\": [\n";
    for (size_t i = 0; i < terms.size(); ++i) {
        const auto& t = terms[i];
        a << "    {\n";
        a << "      \"amplitude\": " << t.amplitude << ",\n";
        a << "      \"frequency\": " << t.frequency << ",\n";
        a << "      \"phase_shift\": " << t.phase_shift << "\n";
        a << "    }" << (i + 1 < terms.size() ? "," : "") << "\n";
    }
    a << "  ],\n\n"
        << "  \"bias\": " << bias << "\n"
        << "}\n";

    // Write merged JSON back
    std::ofstream out(filepath, std::ios::trunc);
    if (!out) {
        std::cerr << "Failed to open solver metadata for writing: " << filepath << "\n";
        return;
    }
    out << content << a.str();
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