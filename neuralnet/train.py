import numpy as np
import joblib

from constants import *
import torch
from nn import cpu_train, gpu_train


if __name__ == "__main__":
    print(f"Loading feature matrices from {FEATURE_MATRICES_PATH}...")
    data = np.load(FEATURE_MATRICES_PATH)
    X, y = data["X"], data["y"]
    print(f"X shape: {X.shape} | y shape: {y.shape}")

    train_fn = gpu_train if torch.cuda.is_available() else cpu_train
    model, X_scaler, y_scaler = train_fn(X, y)

    torch.save(model.state_dict(), NEURAL_NET_DIR / "model.pt")
    joblib.dump(X_scaler, NEURAL_NET_DIR / "X_scaler.pkl")
    joblib.dump(y_scaler, NEURAL_NET_DIR / "y_scaler.pkl")
    print(f"Saved model and scalers to {NEURAL_NET_DIR}")
