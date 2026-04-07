#pragma once

#include <vector>
#include <random>
#include <ostream>
#include <filesystem>

#include "metadata_writer.h"


struct RandomInitialConditionConfig {
    int n = 5;
    double domain_length = 10;
    double amp_min = -1.0;
    double amp_max = 1.0;
    double frequency_multiplier_min = 0.1;
    double frequency_multiplier_max = 5.0;
    int wrap_around_frequency_multiplier_min = 1;
    int wrap_around_frequency_multiplier_max = 5;
};

class RandomInitialCondition : public IMetadataProvider {

public:

    struct Term {
        double amplitude;
        double frequency;
        double phase_shift;
    };

    RandomInitialCondition(
                   const RandomInitialConditionConfig &config,
                   bool alwaysPositive = false,
                   bool wrapAround = true,
                   unsigned seed = std::random_device{}());

    double operator()(double x) const;

    double getMaxU() const;

    void appendMetadata(nlohmann::json& metadata) const override;


    std::string toString() const;

private:
    const RandomInitialConditionConfig& config;
    double domain_length;
    std::vector<Term> terms;
    double bias; // used only if alwaysPositive = true
    double max_u;
};
