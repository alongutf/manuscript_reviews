"""Scratch simulation of a synthetic gene-regulatory network with competitive binding.

Not part of the paper's published pipeline (contrast with src/simulations.py, which
generates the synthetic scRNA-seq data actually used in scripts/simulated_data.ipynb).
This is an exploratory script -- run interactively (e.g. in PyCharm, hence the
Shift+F10 comments below) -- that simulates N_TRAJ independent trajectories of
N_FEATURES "genes" under a random activator/repressor interaction network, each gene
governed by a Michaelis-Menten-like competitive-binding update with linear decay and
additive noise, integrated with a fixed-step Euler scheme.

Reads: nothing.
Writes: results\\sim_test_J-<J>.npy, shape (N_TRAJ, N_FEATURES, N_STEPS), the full
    trajectories in float32. Also plots one gene's trajectory (gene 0) for every
    20th trajectory via plt.show() (interactive; not saved to disk).
Run: python simulations/test.py  (requires a results\\ directory to already exist,
    e.g. run from within simulations/ or wherever the relative path resolves).
"""
# This is a sample Python script.
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
import time
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

@njit
def generate_interaction_matrix(J,n_features):
    """Random gene-gene interaction network, split into an activating and a
    repressing half.

    W's entries are |N(0, (J/sqrt(n_features))^2)|: the 1/sqrt(n_features) scaling
    keeps each gene's total incoming interaction strength O(1) regardless of network
    size, so J alone sets the interaction strength. A is an independent Bernoulli(0.5)
    mask deciding, per entry, whether that interaction acts through the activating
    channel (W_1) or the repressing channel (W_2) in competitive_binding -- so W_1 and
    W_2 partition the same underlying weights rather than being drawn separately.
    """
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
    """Redraw a `ratio` fraction of W1/W2's rows/columns from a fresh network,
    leaving the rest untouched (unused elsewhere in this script; ratio=1 replaces
    everything, equivalent to just calling generate_interaction_matrix again)."""
    k = int(W1.shape[1] * ratio)
    new_W1, new_W2 = generate_interaction_matrix(J, W1.shape[1])
    W1[:k,:] = new_W1[:k,:]
    W2[:k,:] = new_W2[:k,:]
    W1[k:,:k] = new_W1[k:,:k]
    W2[k:,:k] = new_W2[k:,:k]
    return W1, W2


@njit(parallel=True, fastmath=True)
def competitive_binding(n_trajectories, n_features, n_steps, dt, random_state=None):
    """Simulate n_trajectories independent gene-expression time courses.

    Each trajectory starts at x=1 for every gene (arbitrary common initial
    condition) and is integrated forward n_steps with Euler step dt under
        dx/dt = -DEG*x + RATE * (W_act @ x) / (1 + W_rep @ x + W_act @ x) + NOISE * xi,
    i.e. linear degradation, a saturating (competitive-binding-style) activation term
    that is itself damped by the repressive interactions in the denominator, and
    additive Gaussian noise. x is clipped at 0 after every step since expression
    cannot go negative. degradation_rate and noise_amp are read from the enclosing
    module scope (DEG, NOISE) rather than passed in -- numba njit closures over
    module globals here, so changing those globals after this function is compiled
    has no effect until the module is reloaded.

    A fresh random interaction network (W_act, W_rep) is drawn once per trajectory
    (inside the prange loop), so trajectories differ both in the noise realization
    and in the underlying network, not just in noise -- this is what N_TRAJ is
    averaging over downstream when eigenvalue/PCA-style analyses are run on results.

    random_state seeds numpy's global RNG once before generating the noise block;
    because trajectories run in parallel (prange) with no per-trajectory sub-seeding,
    reproducibility with random_state fixed still depends on numba's parallel
    scheduling being deterministic, which is not guaranteed across runs/machines.

    Returns
    -------
    results : float32 array, shape (n_trajectories, n_features, n_steps)
        results[i, :, j] is the expression vector of trajectory i at Euler step j.
    """
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
        # network drawn fresh per trajectory -- the W_act/W_rep computed above
        # (outside the loop) is discarded/unused; see docstring note on N_TRAJ
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

# ---- parameters for this run; see competitive_binding docstring for the model ----
N_TRAJ = 500
N_FEATURES = 1000
N_STEPS = 1000
J= 2                # interaction strength scale (see generate_interaction_matrix)
source=0            # unused
degradation_rate=2  # DEG in competitive_binding (read as a module global, see above)
noise_amp=1         # NOISE in competitive_binding (read as a module global, see above)
start = time.perf_counter()
seeds = [None]       # single un-seeded run; loop scaffolding left in for sweeping seeds
for k,seed in enumerate(seeds):
    results = competitive_binding(N_TRAJ,N_FEATURES, N_STEPS, 0.01, random_state=seed)
    stop = time.perf_counter()
    print(results.shape)
    # plot gene 0's trajectory for every 20th simulated cell/trajectory, as a quick
    # visual sanity check that the dynamics settle rather than blow up or vanish
    for i in range(0,results.shape[0],20):
        for j in range(1):
            plt.plot(results[i,j,:])
    plt.show()
    print(f"Total time: {stop - start:.6f} seconds")
    np.save(f"results\\sim_test_J-{J}.npy", results)