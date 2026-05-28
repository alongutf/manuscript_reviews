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
    # remove diagonal
    W = W - np.diag(W)
    # Split W into two parts based on A
    W_1 = W * A
    W_2 = W * (1 - A)

    return W_1, W_2

@njit
def regenerate_interaction_matrix(W1, W2, J , ratio=1):
    k = int(W1.shape[1] * ratio)
    new_W1, new_W2 = generate_interaction_matrix(J, W1.shape[1])
    W1[:k,:] = new_W1[:k,:]
    W2[:k,:] = new_W2[:k,:]
    W1[k:,:k] = new_W1[k:,:k]
    W2[k:,:k] = new_W2[k:,:k]
    return W1, W2


@njit(parallel=True, fastmath=True)
def competitive_binding(n_trajectories, n_features, n_steps, dt, random_state=None):
    #convert to float32
    dt = np.float32(dt)
    ONE = np.float32(1.0)
    ZERO = np.float32(0.0)
    RATE = np.float32(1.0) # np.random.lognormal(mean=-2,sigma=6,size=n_features).astype(np.float32)
    DEG = np.float32(degradation_rate)
    NOISE = np.float32(noise_amp)
    # 2. Pre-calculate matrices OUTSIDE the loops
    W_act, W_rep = generate_interaction_matrix(J,n_features)

    results = np.zeros((n_trajectories, n_features, n_steps),dtype=np.float32)
    # Pre-allocate random noise block to keep inner loop tight
    # (Optional: doing it inside the loop is also fine in Numba, but this can be cleaner)
    if random_state is not None:
        np.random.seed(random_state)
    noise_block = np.random.randn(n_trajectories, n_steps, n_features).astype(np.float32)

    for i in prange(n_trajectories):
        x = np.ones(n_features, dtype=np.float32)
        W_act, W_rep = generate_interaction_matrix(J, n_features)
        for j in range(n_steps):
            # Calculate the common term v = 1 / (x + 1)
            act_interaction = np.dot(W_act, x)
            rep_interaction = np.dot(W_rep, x)


            # competition term:
            #competition = ONE - x.sum() / CAPACITY

            # Update x
            x += dt * ( -DEG*x + RATE * act_interaction/(ONE+rep_interaction+act_interaction) + NOISE * noise_block[i,j,:])

            # Fast clipping
            x = np.maximum(ZERO, x)

            results[i, :, j] = x

    return results

N_TRAJ = 500
N_FEATURES = 1000
N_STEPS = 1000
J= 2
source=0
degradation_rate=2
noise_amp=1
start = time.perf_counter()
seeds = [None]
for k,seed in enumerate(seeds):
    results = competitive_binding(N_TRAJ,N_FEATURES, N_STEPS, 0.01, random_state=seed)
    stop = time.perf_counter()
    print(results.shape)
    for i in range(0,results.shape[0],20):
        for j in range(1):
            plt.plot(results[i,j,:])
    plt.show()
    print(f"Total time: {stop - start:.6f} seconds")
    np.save(f"results\\sim_test_J-{J}.npy", results)