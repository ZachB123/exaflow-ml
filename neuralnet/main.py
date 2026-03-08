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

from constants import *
from burgers_solution import BurgersSolution


class ViscosityDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32).view(-1, 1)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


class ArtificialViscosityNet(nn.Module):
    def __init__(self, input_dim=5, hidden_dim=128):
        super(ArtificialViscosityNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x):
        return self.net(x)


def get_training_data_folder_names():
    return sorted([item.name for item in DEFAULT_TRAINING_DATA_DIR.iterdir() if item.is_dir() and (item / METADATA_FILENAME).exists()])


def requires_artificial_viscosity(dx, u_i_minus_1, u_i_plus_1):
    ux = (u_i_plus_1 - u_i_minus_1) / (2.0 * dx)
    return ux < 0


# def reverse_engineer_cq(dt, dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu=0, eps=1e-3):
#     uxx = u_i_plus_1 - 2*u_i + u_i_minus_1
#     ux_diff = abs(u_i_plus_1 - u_i_minus_1)
#
#     if abs(uxx) < eps or ux_diff < eps:
#         return None
#
#     return (2 / (dx * ux_diff)) * ((dx**2 / (dt * uxx)) * (u_next_i - u_i + u_i * (dt/dx) * (u_i_plus_1 - u_i_minus_1)) - nu)


def reverse_engineer_cq(dt, dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu=0):
    # todo store nu in the metadata
    # sometimes got divide by 0
    uxx = (u_i_plus_1 - 2 * u_i + u_i_minus_1)
    if uxx == 0:
        return None

    ux_diff = abs(u_i_plus_1 - u_i_minus_1)
    if ux_diff == 0:
        return None

    return (2 / (dx * ux_diff)) * ((dx**2 / (dt * uxx)) * (u_next_i - u_i + u_i * (dt/dx) * (u_i_plus_1 - u_i_minus_1)) - nu)


def get_feature_matrices_for_sample(sample_name):
    X_rows = []
    y_rows = []

    burgers_solution = BurgersSolution(sample_name)
    fine_dt = burgers_solution.time_step_size
    fine_dx = burgers_solution.spatial_step_size
    coarse_dt = fine_dt * COARSENESS_MULTIPLIER
    coarse_dx = fine_dx * COARSENESS_MULTIPLIER
    coarse_num_timesteps = int(burgers_solution.max_time // coarse_dt)
    coarse_num_domain_points = int(burgers_solution.domain_length // coarse_dx)

    # Collect data from multiple timesteps
    for time_step in range(min(coarse_num_timesteps, 50)):
        # just ignoring the boundary condition for now cuz I need to ask yuvi about it
        # and it probably won't make a difference for training
        for spatial_step in range(1, coarse_num_domain_points - 1):
            curr_time = time_step * coarse_dt
            next_time = (time_step + 1) * coarse_dt
            curr_x = spatial_step * coarse_dx
            prev_x = (spatial_step - 1) * coarse_dx
            next_x = (spatial_step + 1) * coarse_dx

            u_i = burgers_solution.get_u(curr_x, curr_time)
            u_next_i = burgers_solution.get_u(curr_x, next_time)
            u_i_minus_1 = burgers_solution.get_u(prev_x, curr_time)
            u_i_plus_1 = burgers_solution.get_u(next_x, curr_time)

            if requires_artificial_viscosity(coarse_dx, u_i_minus_1, u_i_plus_1):
                cq = reverse_engineer_cq(coarse_dt, coarse_dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1)
                if cq is not None and not np.isnan(cq) and not np.isinf(cq):
                    # Clamp cq to reasonable values if necessary
                    cq = np.clip(cq, -100, 100)
                    # we can't use u_next_i as a feature because we will never have that when running the sim normally
                    X_rows.append([coarse_dt, coarse_dx, u_i, u_i_minus_1, u_i_plus_1])
                    y_rows.append(cq)

    X = np.array(X_rows)
    y = np.array(y_rows)
    return X, y


def train_model(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    X_scaler = StandardScaler()
    X_train_scaled = X_scaler.fit_transform(X_train)
    X_test_scaled = X_scaler.transform(X_test)

    train_dataset = ViscosityDataset(X_train_scaled, y_train)
    test_dataset = ViscosityDataset(X_test_scaled, y_test)

    train_loader = DataLoader(train_dataset, batch_size=1024, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=1024, shuffle=False)

    model = ArtificialViscosityNet()
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    epochs = 10
    pbar = tqdm(range(epochs))
    for epoch in pbar:
        model.train()
        train_loss = 0
        for batch_X, batch_y in train_loader:
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        model.eval()
        test_loss = 0
        with torch.no_grad():
            for batch_X, batch_y in test_loader:
                outputs = model(batch_X)
                loss = criterion(outputs, batch_y)
                test_loss += loss.item()

        pbar.set_description(f"Epoch {epoch+1}/{epochs} | Train Loss: {train_loss/len(train_loader):.6f} | Test Loss: {test_loss/len(test_loader):.6f}")

    return model, X_scaler


if __name__ == "__main__":
    samples = get_training_data_folder_names()

    print(f"Found {len(samples)} samples. Extracting features...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(get_feature_matrices_for_sample, samples)

    X_matrices, y_matrices = zip(*results)

    X = np.vstack([m for m in X_matrices if m.size > 0])
    y = np.concatenate([m for m in y_matrices if m.size > 0])

    model, X_scaler = train_model(X, y)

