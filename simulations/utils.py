# Backward-compatibility shim — import directly from src instead of this file.
from src.analysis_functions import scramble, normalize, z_transform, log_transform, get_pcs, get_eig_dist, mp_distribution
from src.data_functions import plot_eig_dist
from src.simulations import calculate_entropy, run_de_analysis, plot_volcano
