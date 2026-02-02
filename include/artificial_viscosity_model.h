#pragma once

#include <torch/torch.h>
#include <vector>

class ArtificialViscosityModel : public torch::nn::Module {
public:
    // Constructor that takes a vector of hidden layer sizes
    // e.g. [32, 32, 16] creates input(5) -> 32 -> 32 -> 16 -> output(1)
    ArtificialViscosityModel(const std::vector<int>& hidden_layer_sizes);

    // Forward pass: takes 5 inputs and returns 1 output
    // Input: [u_i-1, u_i, u_i+1, first_derivative, second_derivative]
    // Output: predicted artificial viscosity term at position i
    torch::Tensor forward(torch::Tensor x);

private:
    torch::nn::Sequential layers;
};