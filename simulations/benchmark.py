import numpy as np
from numba import njit, prange
import time

J = 0.5


# --- VERSION 1: FLOAT 64 (Standard) ---
@njit(parallel=True, fastmath=True)
def run_f64(n_trajectories, n_features, n_steps, dt):
    # Standard Python floats are 64-bit
    dt = np.float64(dt)

    # Pre-calculate matrices (Float64)
    W = (J / np.sqrt(n_features)) * np.random.randn(n_features, n_features)
    A = np.random.randint(0, 2, size=(n_features, n_features)).astype(np.float64)

    W_1 = W * A
    W_2 = W * (1 - A)
    M_diff = W_2 - W_1
    bias_vec = np.sum(W_1, axis=1)

    results = np.zeros((n_trajectories, n_steps, n_features), dtype=np.float64)

    for i in prange(n_trajectories):
        x = np.ones(n_features, dtype=np.float64)
        noise_block = np.random.randn(n_steps, n_features)  # Default is float64

        for j in range(n_steps):
            v = 1.0 / (x + 1.0)
            interaction = bias_vec + (M_diff @ v)
            x += dt * (-x + interaction + noise_block[j])
            x = np.maximum(0.0, x)
            results[i, j, :] = x

    return results


# --- VERSION 2: FLOAT 32 (Optimized) ---
@njit(parallel=True, fastmath=True)
def run_f32(n_trajectories, n_features, n_steps, dt):
    # 1. Define explicit 32-bit constants
    dt = np.float32(dt)
    ONE = np.float32(1.0)
    ZERO = np.float32(0.0)

    # Cast matrices to float32
    W = (J / np.sqrt(n_features)) * np.random.randn(n_features, n_features).astype(np.float32)
    A = np.random.randint(0, 2, size=(n_features, n_features)).astype(np.float32)

    W_1 = W * A
    W_2 = W * (1 - A)
    M_diff = (W_2 - W_1).astype(np.float32)
    bias_vec = np.sum(W_1, axis=1).astype(np.float32)

    results = np.zeros((n_trajectories, n_steps, n_features), dtype=np.float32)

    for i in prange(n_trajectories):
        x = np.ones(n_features, dtype=np.float32)
        noise_block = np.random.randn(n_steps, n_features).astype(np.float32)

        for j in range(n_steps):
            # ERROR WAS HERE: "1.0" promoted this to float64
            # FIX: Use "ONE" (which is float32)
            v = ONE / (x + ONE)

            interaction = bias_vec + (M_diff @ v)
            x += dt * (-x + interaction + noise_block[j])

            # ERROR WAS HERE: "0.0" promoted this to float64
            # FIX: Use "ZERO" (which is float32)
            x = np.maximum(ZERO, x)

            results[i, j, :] = x

    return results

# --- BENCHMARK HARNESS ---
def benchmark():
    # Parameters (Scaled down slightly for quick testing)
    N_TRAJ = 100
    N_FEAT = 1000
    N_STEP = 1000
    DT = 0.01

    print(f"Benchmarking: {N_TRAJ} trajectories, {N_FEAT} features, {N_STEP} steps")
    print("-" * 50)

    # 1. Warmup (Compile)
    print("Compiling Float64...")
    _ = run_f64(10, 100, 10, 0.01)
    print("Compiling Float32...")
    _ = run_f32(10, 100, 10, 0.01)

    # 2. Run Float64
    print("Running Float64 simulation...")
    start = time.perf_counter()
    _ = run_f64(N_TRAJ, N_FEAT, N_STEP, DT)
    end = time.perf_counter()
    t64 = end - start
    print(f"Float64 Time: {t64:.4f} seconds")

    # 3. Run Float32
    print("Running Float32 simulation...")
    start = time.perf_counter()
    _ = run_f32(N_TRAJ, N_FEAT, N_STEP, DT)
    end = time.perf_counter()
    t32 = end - start
    print(f"Float32 Time: {t32:.4f} seconds")

    # 4. Results
    print("-" * 50)
    print(f"Speedup: {t64 / t32:.2f}x faster")


if __name__ == "__main__":
    benchmark()