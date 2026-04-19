# import torch
# from nn import ArtificialViscosityNet
# from constants import U_RADIUS

# # Reconstruct model and load existing weights
# model = ArtificialViscosityNet()
# model.load_state_dict(torch.load("output/model.pt", map_location=torch.device("cpu")))
# model.eval()

# # Trace and save as TorchScript
# example_input = torch.randn(1, 6 + U_RADIUS * 2)
# scripted = torch.jit.trace(model, example_input)
# scripted.save("output/model_scripted.pt")
# print("Done")

import torch
import joblib
import numpy as np
from nn import ArtificialViscosityNet
from constants import U_RADIUS

model = torch.jit.load("output/model.pt")
model.eval()

X_scaler = joblib.load("output/X_scaler.pkl")
y_scaler = joblib.load("output/y_scaler.pkl")

# Bake scaler params into the model as a wrapper
class ScaledModel(torch.nn.Module):
    def __init__(self, model, X_scaler, y_scaler):
        super().__init__()
        self.model = model
        self.register_buffer("X_mean", torch.tensor(X_scaler.mean_, dtype=torch.float32))
        self.register_buffer("X_std",  torch.tensor(X_scaler.scale_, dtype=torch.float32))
        self.register_buffer("y_mean", torch.tensor(y_scaler.mean_, dtype=torch.float32))
        self.register_buffer("y_std",  torch.tensor(y_scaler.scale_, dtype=torch.float32))

    def forward(self, x):
        x = (x - self.X_mean) / self.X_std
        y = self.model(x)
        return y * self.y_std + self.y_mean

scaled = ScaledModel(model, X_scaler, y_scaler)
scaled.eval()

example = torch.randn(1, 6 + U_RADIUS * 2)
scripted = torch.jit.trace(scaled, example)
scripted.save("output/model.pt")
print("Done")