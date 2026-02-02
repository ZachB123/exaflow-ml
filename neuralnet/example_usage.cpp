#include "artificial_viscosity_model.h"
#include <torch/torch.h>
#include <vector>

// DON'T ACTUALLY RUN THIS, IT IS JUST AN EXAMPLE

// Create model
std::vector<int> hidden_sizes = {32, 32, 16};
auto model = std::make_shared<ArtificialViscosityModel>(hidden_sizes);

// Training loop example
torch::optim::Adam optimizer(model->parameters(), 0.001);
for (auto& batch : training_data) {
    optimizer.zero_grad();
    auto output = model->forward(batch.input);
    auto loss = torch::nn::functional::mse_loss(output, batch.target);
    loss.backward();
    optimizer.step();
}

// Save/load model weights
torch::save(model, "model_weights.pt");
torch::load(model, "model_weights.pt");

// Prediction
model->eval();  // Set to evaluation mode
torch::NoGradGuard no_grad;  // Disable gradient computation
// Input: [u_i-1, u_i, u_i+1, first_derivative, second_derivative]
auto input_tensor = torch::tensor({{1.0, 2.0, 1.5, -0.5, 0.1}});
auto prediction = model->forward(input_tensor);