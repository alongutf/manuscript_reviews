# This is a sample Python script.
"""
Standalone gene-expression-dynamics simulator (competition-free, numba-accelerated).

Simulates N_TRAJ independent trajectories of an interaction network of N_FEATURES
"genes", each governed by a linear degradation term, a constant production rate, a
pairwise interaction term through a random signed weight matrix, and additive noise.
This is not one of the GMP-Cor / eigenvalue-spectrum runners described in CLAUDE.md's
"simulations/" table (it predates that convention and does not write to
results/simulation_results/); it is a dynamical-systems exploration script that saves
raw trajectory arrays for offline inspection.

Reads: nothing (all parameters are hard-coded in the __main__ block below).
Writes: results\\sim4_J-{J}_rep{k+1}.npy — one array per entry of `seeds`, shape
    (n_trajectories, n_features, n_steps), float32. The path is relative to the
    current working directory, not repo-root-anchored.
Run: python simulations/main.py  (from the repo root, so the relative "results\\..."
    path resolves to the existing top-level results/ directory — see FINDINGS in the
    accompanying log for what happens if that assumption doesn't hold).
"""
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
import time
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

@njit
def generate_interaction_matrix(J,n_features):
    # Build the random gene-gene interaction matrix and its constant bias term.
    #
    # J sets the interaction scale: entries of W are |N(0, (J/sqrt(n_features))^2)|,
    # so all interactions are non-negative before the sign is assigned by A (below).
    # The sqrt(n_features) scaling keeps the total interaction strength per gene
    # roughly independent of network size (standard random-matrix normalization).
    # A is an independent 0/1 coin per entry that decides, per gene pair, whether that
    # interaction acts through the "W_2 - v" branch or the "W_1 - (1-v)" branch of the
    # underlying rate equation (see the simplification note below); it is NOT a sign of
    # W itself.
    #
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
    # bias = W_1 @ 1 (a vector of ones), i.e. the row-sums of W_1 -- the "(1-v)=1" term
    # of the original two-matmul formula, meant to be added back on every step of the
    # simulation loop (see run_simulation_optimized's `interaction` line).
    return M_diff, np.sum(W_1, axis=1).astype(np.float32)


@njit(parallel=True, fastmath=True)
def run_simulation_optimized(n_trajectories, n_features, n_steps, dt, weights=None,bias=None, random_state=None):
    # Integrate n_trajectories independent copies of the gene network in parallel
    # (prange splits trajectories across threads), each starting from x = 1 for every
    # gene, using a fixed-step explicit Euler update:
    #   dx/dt = -DEG*x + RATE + interaction(x) + NOISE * dW
    # with x clipped at 0 (concentrations/counts cannot go negative) after every step.
    #
    # degradation_rate, noise_amp and J are read here as module-level globals (not
    # function parameters) -- they must already be set at module scope before this
    # jitted function is first called and traced/compiled, and numba's nopython mode
    # treats such globals as compile-time constants: reassigning them after the first
    # call will NOT change the compiled function's behaviour.
    #convert to float32
    dt = np.float32(dt)
    ONE = np.float32(1.0)
    ZERO = np.float32(0.0)
    RATE = np.float32(10.0) # np.random.lognormal(mean=-2,sigma=6,size=n_features).astype(np.float32)
    DEG = np.float32(degradation_rate)
    NOISE = np.float32(noise_amp)
    CAPACITY = np.float32(5*n_features)  # unused: the logistic competition term below is disabled
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

            # Optimized update step: 1 matrix mult instead of 2.
            # Per the algebraic-simplification note in generate_interaction_matrix,
            # the full simplified interaction term is (weights @ v) + bias; `bias` is
            # accepted as a parameter and returned by generate_interaction_matrix but is
            # never added here -- see the log file's FINDINGS section.
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
    N_TRAJ = 500       # number of independent trajectories to simulate
    N_FEATURES = 1000  # number of interacting "genes" per trajectory
    N_STEPS = 1000      # number of Euler integration steps
    J=5                 # interaction-strength scale, see generate_interaction_matrix
    source=0            # unused
    degradation_rate=1  # DEG: linear decay rate, read as a global inside the njit function
    noise_amp=3         # NOISE: additive noise amplitude, read as a global inside the njit function
    start = time.perf_counter()
    interaction_matrix, bias_vec = generate_interaction_matrix(J,N_FEATURES)
    print(np.mean(bias_vec))
    seeds = [None]  # single run, no fixed RNG seed -> not reproducible run-to-run
    for k,seed in enumerate(seeds):
        results = run_simulation_optimized(N_TRAJ,N_FEATURES, N_STEPS, 0.01, weights=interaction_matrix, bias=bias_vec, random_state=seed)
        stop = time.perf_counter()
        print(results.shape)
        # Overlay gene 0's trajectory for every 100th simulated trajectory
        for i in range(0,results.shape[0],100):
            for j in range(1):
                plt.plot(results[i,j,:])
        plt.show()
        print(f"Total time: {stop - start:.6f} seconds")
        # Path is relative to the current working directory (not repo-root-anchored) and
        # assumes a "results" directory already exists there -- see FINDINGS in the log.
        np.save(f"results\\sim4_J-{J}_rep{k+1}.npy", results)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
