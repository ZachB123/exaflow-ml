import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from tqdm import tqdm
from torch.utils.data import DataLoader, Dataset
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import joblib
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

from constants import *

def get_device():
    if torch.cuda.is_available():
        device = torch.device("cuda")
    # elif torch.backends.mps.is_available():
    #     device = torch.device("mps")
    else:
        device = torch.device("cpu")
    print(f"Using device: {device}")
    return device

def plot_losses(train_losses, test_losses):
    plt.figure(figsize=(8, 5))
    plt.plot(train_losses, label="Train Loss")
    plt.plot(test_losses, label="Test Loss")
    plt.xlabel("Epoch")
    plt.ylabel("MSE Loss")
    plt.yscale("log")
    plt.legend()
    plt.title("Training & Test Loss")
    plt.tight_layout()
    plt.savefig("losses.png", dpi=150)
    plt.show()


class ViscosityDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32).view(-1, 1)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


class ArtificialViscosityNet(nn.Module):
    # 6 features for dx, dt, ux, uxx, u_curr and nu, then double the radius features
    def __init__(self, input_dim=(6 + U_RADIUS * 2), hidden_dim=256):
        super(ArtificialViscosityNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x):
        return self.net(x)


def cpu_train(X, y):
    device = get_device()

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    X_scaler = StandardScaler()
    X_train_scaled = X_scaler.fit_transform(X_train)
    X_test_scaled = X_scaler.transform(X_test)

    y_scaler = StandardScaler()
    y_train_scaled = y_scaler.fit_transform(y_train.reshape(-1, 1)).ravel()
    y_test_scaled = y_scaler.transform(y_test.reshape(-1, 1)).ravel()

    train_dataset = ViscosityDataset(X_train_scaled, y_train_scaled)
    test_dataset = ViscosityDataset(X_test_scaled, y_test_scaled)

    # pin_memory=False: no benefit on CPU (only useful when sending data to a GPU)
    # num_workers=cpu_count(): let all cores share the data loading work
    train_loader = DataLoader(train_dataset, batch_size=2048, shuffle=True,
                              num_workers=cpu_count(), pin_memory=False)
    test_loader = DataLoader(test_dataset, batch_size=2048, shuffle=False,
                             num_workers=cpu_count(), pin_memory=False)

    model = ArtificialViscosityNet().to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10)

    train_losses, test_losses = [], []
    pbar = tqdm(range(EPOCHS))
    for epoch in pbar:
        model.train()
        train_loss = 0
        for batch_X, batch_y in train_loader:
            optimizer.zero_grad()
            loss = criterion(model(batch_X), batch_y)
            loss.backward()
            optimizer.step()
            train_loss += loss

        model.eval()
        test_loss = 0
        with torch.no_grad():
            for batch_X, batch_y in test_loader:
                test_loss += criterion(model(batch_X), batch_y)

        avg_train = (train_loss / len(train_loader)).item()
        avg_test = (test_loss / len(test_loader)).item()
        train_losses.append(avg_train)
        test_losses.append(avg_test)
        scheduler.step(avg_test)
        pbar.set_description(f"Epoch {epoch+1}/{EPOCHS} | Train: {avg_train:.6f} | Test: {avg_test:.6f}")

    print(train_losses)
    print()
    print(test_losses)
    plot_losses(train_losses, test_losses)

    return model, X_scaler, y_scaler


def gpu_train(X, y):
    device = get_device()
    assert device.type == "cuda", "gpu_train requires a CUDA device"

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    X_scaler = StandardScaler()
    X_train_scaled = X_scaler.fit_transform(X_train)
    X_test_scaled = X_scaler.transform(X_test)

    y_scaler = StandardScaler()
    y_train_scaled = y_scaler.fit_transform(y_train.reshape(-1, 1)).ravel()
    y_test_scaled = y_scaler.transform(y_test.reshape(-1, 1)).ravel()

    # Move full dataset to GPU once — eliminates all CPU→GPU transfers during training
    X_train_t = torch.tensor(X_train_scaled, dtype=torch.float32).to(device)
    y_train_t = torch.tensor(y_train_scaled, dtype=torch.float32).view(-1, 1).to(device)
    X_test_t  = torch.tensor(X_test_scaled,  dtype=torch.float32).to(device)
    y_test_t  = torch.tensor(y_test_scaled,  dtype=torch.float32).view(-1, 1).to(device)
    n_train, n_test = X_train_t.shape[0], X_test_t.shape[0]

    BATCH_SIZE = 65536

    model = ArtificialViscosityNet().to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10)

    train_losses, test_losses = [], []
    pbar = tqdm(range(EPOCHS))
    for epoch in pbar:
        model.train()
        train_loss = 0
        n_train_batches = 0
        perm = torch.randperm(n_train, device=device)
        for start in range(0, n_train, BATCH_SIZE):
            idx = perm[start:start + BATCH_SIZE]
            optimizer.zero_grad()
            loss = criterion(model(X_train_t[idx]), y_train_t[idx])
            loss.backward()
            optimizer.step()
            train_loss += loss
            n_train_batches += 1

        model.eval()
        test_loss = 0
        n_test_batches = 0
        with torch.no_grad():
            for start in range(0, n_test, BATCH_SIZE):
                test_loss += criterion(model(X_test_t[start:start + BATCH_SIZE]),
                                       y_test_t[start:start + BATCH_SIZE])
                n_test_batches += 1

        avg_train = (train_loss / n_train_batches).item()
        avg_test = (test_loss / n_test_batches).item()
        train_losses.append(avg_train)
        test_losses.append(avg_test)
        scheduler.step(avg_test)
        pbar.set_description(f"Epoch {epoch+1}/{EPOCHS} | Train: {avg_train:.6f} | Test: {avg_test:.6f}")

    print(train_losses)
    print()
    print(test_losses)
    plot_losses(train_losses, test_losses)

    return model, X_scaler, y_scaler


