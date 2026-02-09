#include <iostream>
#include <random>
#include <filesystem>
#include <memory>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>

#include "initial_condition_generator.h"
#include "burgers.h"

// how many training samples to generate
const int DEFAULT_NUM_SAMPLES = 2;
const std::string TRAINING_DIR = "../training_data";

// what stencil we are using for the configuration
const auto STENCIL_FACTORY = []() {
    return std::make_unique<LaxWendroff>();
};

// constant solver configurations
const double KINEMATIC_VISCOSITY = 0.01;
const int DEFAULT_TIME_STEPS = 10000;
const double DEFAULT_TIME_STEP_SIZE = 0.000001;

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

int findNextSampleNumber() {
    int max_sample = -1;
    
    if (!std::filesystem::exists(TRAINING_DIR)) {
        return 0;
    }
    
    for (const auto& entry : std::filesystem::directory_iterator(TRAINING_DIR)) {
        if (entry.is_directory()) {
            std::string folder_name = entry.path().filename().string();
            if (folder_name.substr(0, 7) == "sample_") {
                try {
                    int sample_num = std::stoi(folder_name.substr(7));
                    max_sample = std::max(max_sample, sample_num);
                } catch (...) {
                    // ignore folders that don't match the pattern
                }
            }
        }
    }
    
    return max_sample + 1;
}

void deleteExistingSamples() {
    if (std::filesystem::exists(TRAINING_DIR)) {
        std::cout << "Clearing existing training data...\n";
        for (const auto& entry : std::filesystem::directory_iterator(TRAINING_DIR)) {
            if (entry.is_directory() && entry.path().filename().string().substr(0, 7) == "sample_") {
                std::filesystem::remove_all(entry.path());
            }
        }
    }
}

void parseCommandLineArguments(int argc, char* argv[], int& num_samples, bool& append_mode, int& time_steps, double& time_step_size) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--samples" || arg == "-s") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << arg << " requires a value\n";
                std::exit(1);
            }
            try {
                num_samples = std::stoi(argv[++i]);
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid value for " << arg << ": " << argv[i] << "\n";
                std::exit(1);
            }
        }
        else if (arg == "--append" || arg == "-a") {
            append_mode = true;
        }
        else if (arg == "--time-step" || arg == "-ts") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << arg << " requires a value\n";
                std::exit(1);
            }
            try {
                time_step_size = std::stod(argv[++i]);
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid value for " << arg << ": " << argv[i] << "\n";
                std::exit(1);
            }
        }
        else if (arg == "--num-time-steps" || arg == "-nts") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << arg << " requires a value\n";
                std::exit(1);
            }
            try {
                time_steps = std::stoi(argv[++i]);
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid value for " << arg << ": " << argv[i] << "\n";
                std::exit(1);
            }
        }
        else {
            std::cerr << "Error: Unknown argument '" << arg << "'\n";
            std::cerr << "Valid arguments:\n";
            std::cerr << "  --samples, -s <num>          Number of samples to generate\n";
            std::cerr << "  --append, -a                 Append to existing samples\n";
            std::cerr << "  --time-step, -ts <value>     Time step size\n";
            std::cerr << "  --num-time-steps, -nts <num> Number of time steps\n";
            std::exit(1);
        }
    }
}

int main(int argc, char* argv[]) {

    int num_samples = DEFAULT_NUM_SAMPLES;
    bool append_mode = false;
    int time_steps = DEFAULT_TIME_STEPS;
    double time_step_size = DEFAULT_TIME_STEP_SIZE;
    
    parseCommandLineArguments(argc, argv, num_samples, append_mode, time_steps, time_step_size);
    
    int start_sample_index = 0;
    
    if (append_mode) {
        start_sample_index = findNextSampleNumber();
        std::cout << "Append mode: starting from sample " << start_sample_index << "\n";
    } else {
        deleteExistingSamples();
        std::cout << "Starting from sample 0\n";
    }
    
    std::cout << "Generating " << num_samples << " samples with " 
              << time_steps << " time steps of size " << time_step_size << "\n";

    for (int sample_index = start_sample_index; sample_index < start_sample_index + num_samples; ++sample_index) {

        RandomInitialConditionConfig cfg = generateRandomInitialConditionConfig();

        RandomInitialCondition f(cfg);

        SolverConfig solver_cfg;
        solver_cfg.kinematic_viscosity = KINEMATIC_VISCOSITY;
        solver_cfg.time_steps = time_steps;
        solver_cfg.domain_length = cfg.domain_length;
        solver_cfg.time_step_size = time_step_size;

        BurgersSolver1d solver(
            STENCIL_FACTORY(),
            solver_cfg,
            f
        );

        solver.solve();
        std::cout << "random function was nan detected: " << solver.wasNanDetected() << std::endl;

        std::ostringstream folder_name;
        folder_name << "sample_" << std::setw(6) << std::setfill('0') << sample_index;

        solver.saveSolution(TRAINING_DIR, folder_name.str(), 1);
        f.saveMetadataJSON(TRAINING_DIR, folder_name.str(),
                          time_steps, time_step_size,
                          solver.getNumDomainPoints(),
                          solver.getSpatialStepSize(),
                          solver.getStencilName());

        std::cout << "Generated sample " << sample_index << " in "
                  << TRAINING_DIR << "/" << folder_name.str() << "\n";
    }

    return 0;
}
