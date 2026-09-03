# Standalone exploratory simulation of a "competitive binding" gene circuit: each
# gene's rate of change depends on a weighted sum of activating and repressing
# inputs from other genes, run forward as a stochastic ODE with numba-jitted inner
# loops for speed. This is not wired into the AnnMat/GMP-Cor pipeline described in
# CLAUDE.md and is not referenced by any other script in the repo (leftover
# IDE-template comment below confirms it was a scratch file, not a maintained one).
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
import time
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

@njit
def generate_interaction_matrix(J,n_features):
    # We generate standard normal noise once
    # |N(0, (J/sqrt(n_features))^2)|: a non-negative random weight matrix, scaled by
    # 1/sqrt(n_features) so total input to a gene stays O(1) as n_features grows
    W = (np.abs((J / np.sqrt(n_features)) * np.random.randn(n_features, n_features))).astype(np.float32)
    # A is a random 0/1 mask that splits every edge into "activating" (W_1) or
    # "repressing" (W_2); the same edge cannot be both, so A partitions W exactly
    A = np.random.randint(0, 2, size=(n_features, n_features)).astype(np.float32)
    # remove diagonal
    W = W - np.diag(W)
    # Split W into two parts based on A
    W_1 = W * A
    W_2 = W * (1 - A)

    return W_1, W_2

@njit
def regenerate_interaction_matrix(W1, W2, J , ratio=1):
    # Rewires a fraction `ratio` of the first k rows (and the corresponding block of
    # the remaining rows' first k columns) of an existing interaction matrix, drawing
    # fresh weights via generate_interaction_matrix and splicing them in in place.
    # Not called anywhere in this script (defined but unused) — a leftover helper
    # for an experiment that partially re-randomizes the network between runs.
    k = int(W1.shape[1] * ratio)
    new_W1, new_W2 = generate_interaction_matrix(J, W1.shape[1])
    W1[:k,:] = new_W1[:k,:]
    W2[:k,:] = new_W2[:k,:]
    W1[k:,:k] = new_W1[k:,:k]
    W2[k:,:k] = new_W2[k:,:k]
    return W1, W2


@njit(parallel=True, fastmath=True)
def competitive_binding(n_trajectories, n_features, n_steps, dt, random_state=None):
    # Simulates n_trajectories independent gene-expression trajectories under a
    # Langevin-type update: deterministic decay + a saturating activation/repression
    # term (activators and repressors compete for the same denominator, hence the
    # name) + additive Gaussian noise. Reads J, degradation_rate and noise_amp from
    # module-level globals (set below, before this function is first called) rather
    # than from its own arguments — numba bakes those in as compile-time constants
    # at first call, so changing the globals after the first call has no effect.
    # Returns an array of shape (n_trajectories, n_features, n_steps): one gene
    # expression trace per trajectory per gene per time step.
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
        # each trajectory gets its own random interaction network, not the one
        # precomputed above — see FINDINGS in the review log for the consequence
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

# ── Parameters and run ──────────────────────────────────────────────────────
N_TRAJ = 500
N_FEATURES = 1000
N_STEPS = 1000
J= 2.5                  # interaction strength scale fed into generate_interaction_matrix
source=0                # unused — left over from an earlier version of the model
degradation_rate=2      # per-step decay rate DEG in competitive_binding
noise_amp=1             # additive Gaussian noise amplitude NOISE in competitive_binding
start = time.perf_counter()
seeds = [None]           # single unseeded run; extend this list to compare seeds
for k,seed in enumerate(seeds):
    results = competitive_binding(N_TRAJ,N_FEATURES, N_STEPS, 0.01, random_state=seed)
    stop = time.perf_counter()
    print(results.shape)
    # plot gene 0's trajectory for every 20th simulated trajectory, as a quick
    # visual sanity check that the dynamics settle rather than diverge or vanish
    for i in range(0,results.shape[0],20):
        for j in range(1):
            plt.plot(results[i,j,:])
    plt.show()
    print(f"Total time: {stop - start:.6f} seconds")
    # relative path: only correct when this script is run with the repo root as
    # the working directory (a results\ folder must already exist there)
    np.save(f"results\\sim5_competitive_J-{J}_regenerate.npy", results)