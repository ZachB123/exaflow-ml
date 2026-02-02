#include <iostream>
#include <random>
#include <filesystem>
#include <memory>
#include <iomanip>
#include <sstream>

#include "initial_condition_generator.h"
#include "burgers.h"

// how many training samples to generate
const int NUM_SAMPLES = 5;
const std::string TRAINING_DIR = "../training_data"; 

// constant solver configurations
const double KINEMATIC_VISCOSITY = 0.01;
const int TIME_STEPS = 1000; // small for debugging purposes
const double TIME_STEP_SIZE = 0.005; // large for debugging purposes

// n (number of terms)
const int N_MIN = 1;
const int N_MAX = 20;

// domain length
const double DOMAIN_MIN = 2.0;
const double DOMAIN_MAX = 100.0;

// amplitude range
const double AMP_MIN = -10.0;
const double AMP_MAX = 10.0;

// frequency multiplier range
const double FREQ_MIN = 0.1;
const double FREQ_MAX = 20.0;

// wrap-around frequency multiplier
const int WRAP_K_MIN = 1;
const int WRAP_K_MAX = 20;


RandomInitialConditionConfig generateRandomInitialConditionConfig() {
    static std::random_device rd;
    static std::mt19937 rng(rd());

    std::uniform_int_distribution<int> n_dist(N_MIN, N_MAX);
    std::uniform_real_distribution<double> domain_dist(DOMAIN_MIN, DOMAIN_MAX);
    std::uniform_real_distribution<double> amp_dist(AMP_MIN, AMP_MAX);
    std::uniform_real_distribution<double> freq_dist(FREQ_MIN, FREQ_MAX);
    std::uniform_int_distribution<int> wrap_k_dist(WRAP_K_MIN, WRAP_K_MAX);

    RandomInitialConditionConfig cfg;

    cfg.n = n_dist(rng);
    cfg.domain_length = domain_dist(rng);

    // Amplitude
    double a1 = amp_dist(rng);
    double a2 = amp_dist(rng);
    cfg.amp_min = std::min(a1, a2);
    cfg.amp_max = std::max(a1, a2);

    // Frequency multiplier
    double f1 = freq_dist(rng);
    double f2 = freq_dist(rng);
    cfg.frequency_multiplier_min = std::min(f1, f2);
    cfg.frequency_multiplier_max = std::max(f1, f2);

    // Wrap-around frequency multiplier
    int w1 = wrap_k_dist(rng);
    int w2 = wrap_k_dist(rng);
    cfg.wrap_around_frequency_multiplier_min = std::min(w1, w2);
    cfg.wrap_around_frequency_multiplier_max = std::max(w1, w2);

    return cfg;
}


int main() {

    for (int sample_index = 0; sample_index < NUM_SAMPLES; ++sample_index) {

        RandomInitialConditionConfig cfg = generateRandomInitialConditionConfig();

        RandomInitialCondition f(cfg);

        SolverConfig solver_cfg;
        solver_cfg.kinematic_viscosity = KINEMATIC_VISCOSITY;
        solver_cfg.time_steps = TIME_STEPS;
        solver_cfg.domain_length = cfg.domain_length;
        solver_cfg.time_step_size = TIME_STEP_SIZE;

        BurgersSolver1d solver(
            std::make_unique<LaxWendroff>(),
            solver_cfg,
            f
        );

        solver.solve();
        std::cout << "random function was nan detected: " << solver.wasNanDetected() << std::endl;

        std::ostringstream folder_name;
        folder_name << "sample_" << std::setw(6) << std::setfill('0') << sample_index;

        solver.saveSolution(TRAINING_DIR, folder_name.str(), 1);
        f.saveMetadataJSON(TRAINING_DIR, folder_name.str());

        std::cout << "Generated sample " << sample_index << " in "
                  << TRAINING_DIR << "/" << folder_name.str() << "\n";
    }

    return 0;
}
