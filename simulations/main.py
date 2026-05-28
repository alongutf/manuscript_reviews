# This is a sample Python script.
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
import time
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

@njit
def generate_interaction_matrix(J,n_features):
    # We generate standard normal noise once
    W = (np.abs((J / np.sqrt(n_features)) * np.random.randn(n_features, n_features))).astype(np.float32)
    A = np.random.randint(0, 2, size=(n_features, n_features)).astype(np.float32)

    # Split W into two parts based on A
    W_1 = W * A
    W_2 = W * (1 - A)

    # 3. Algebraic Simplification
    # Original logic: W_1 @ (1-v) + W_2 @ v
    # Simplified: (W_2 - W_1) @ v + (W_1 @ 1)
    # This reduces 2 matrix-mults to 1 matrix-mult per step.
    M_diff = (W_2 - W_1).astype(np.float32)  # The combined matrix
    return M_diff, np.sum(W_1, axis=1).astype(np.float32)


@njit(parallel=True, fastmath=True)
def run_simulation_optimized(n_trajectories, n_features, n_steps, dt, weights=None,bias=None, random_state=None):
    #convert to float32
    dt = np.float32(dt)
    ONE = np.float32(1.0)
    ZERO = np.float32(0.0)
    RATE = np.float32(10.0) # np.random.lognormal(mean=-2,sigma=6,size=n_features).astype(np.float32)
    DEG = np.float32(degradation_rate)
    NOISE = np.float32(noise_amp)
    CAPACITY = np.float32(5*n_features)
    # 2. Pre-calculate matrices OUTSIDE the loops
    if weights is None:
        weights, bias = generate_interaction_matrix(J,n_features)

    results = np.zeros((n_trajectories, n_features, n_steps),dtype=np.float32)
    # Pre-allocate random noise block to keep inner loop tight
    # (Optional: doing it inside the loop is also fine in Numba, but this can be cleaner)
    if random_state is not None:
        np.random.seed(random_state)
    noise_block = np.random.randn(n_trajectories, n_steps, n_features).astype(np.float32)

    for i in prange(n_trajectories):
        x = np.ones(n_features, dtype=np.float32)

        for j in range(n_steps):
            # Calculate the common term v = 1 / (x + 1)
            v = ONE / (x + ONE)

            # Optimized update step: 1 matrix mult instead of 2
            interaction = weights @ v

            # competition term:
            #competition = ONE - x.sum() / CAPACITY

            # Update x
            x += dt * ( -DEG*x + RATE + interaction + NOISE * noise_block[i,j,:])

            # Fast clipping
            x = np.maximum(ZERO, x)

            results[i, :, j] = x

    return results


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    N_TRAJ = 500
    N_FEATURES = 1000
    N_STEPS = 1000
    J=5
    source=0
    degradation_rate=1
    noise_amp=3
    start = time.perf_counter()
    interaction_matrix, bias_vec = generate_interaction_matrix(J,N_FEATURES)
    print(np.mean(bias_vec))
    seeds = [None]
    for k,seed in enumerate(seeds):
        results = run_simulation_optimized(N_TRAJ,N_FEATURES, N_STEPS, 0.01, weights=interaction_matrix, bias=bias_vec, random_state=seed)
        stop = time.perf_counter()
        print(results.shape)
        for i in range(0,results.shape[0],100):
            for j in range(1):
                plt.plot(results[i,j,:])
        plt.show()
        print(f"Total time: {stop - start:.6f} seconds")
        np.save(f"results\\sim4_J-{J}_rep{k+1}.npy", results)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
