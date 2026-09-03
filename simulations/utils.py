# Backward-compatibility shim — import directly from src instead of this file.
# Re-exports the eigenvalue/GMP-Cor and simulation helpers under simulations.utils
# so old scripts/notebooks written against this path keep working unmodified;
# nothing here is implemented locally, it all lives in the src/ modules.
from src.analysis_functions import scramble, normalize, z_transform, log_transform, get_pcs, get_eig_dist, mp_distribution
from src.data_functions import plot_eig_dist
from src.simulations import calculate_entropy, run_de_analysis, plot_volcano
