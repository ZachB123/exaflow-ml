#include "artificial_viscosity_model.h"

ArtificialViscosityModel::ArtificialViscosityModel(const std::vector<int>& hidden_layer_sizes) {
    const int input_size = 5;  // u_i-1, u_i, u_i+1, first_derivative, second_derivative
    const int output_size = 1; // predicted artificial viscosity term

    // Build the sequential layers
    torch::nn::Sequential seq;

    // Connecting the sequential layers
    if (hidden_layer_sizes.empty()) {
        // If no hidden layers specified, create a simple 5 -> 1 network
        seq->push_back(torch::nn::Linear(input_size, output_size));
    } else {
        // Input to first hidden layer
        seq->push_back(torch::nn::Linear(input_size, hidden_layer_sizes[0]));
        seq->push_back(torch::nn::ReLU());

        // Hidden layers with ReLU activation
        for (size_t i = 1; i < hidden_layer_sizes.size(); ++i) {
            seq->push_back(torch::nn::Linear(hidden_layer_sizes[i - 1], hidden_layer_sizes[i]));
            seq->push_back(torch::nn::ReLU());
        }

        // Final hidden layer to output
        seq->push_back(torch::nn::Linear(hidden_layer_sizes.back(), output_size));
    }

    // Register the sequential module
    layers = register_module("layers", seq);
}

torch::Tensor ArtificialViscosityModel::forward(torch::Tensor x) {
    return layers->forward(x);
}