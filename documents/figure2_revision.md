# revision plan for figure2
figure 2 is the main conceptual introduction to the idea of understanding global correlation strength from correlation eigenvalues.

**main points that need to be revised**:

- instead of putting the focus on the GMP analytical distribution, I need to move the spotlight to the synthetic data simulation.
- a schematic graphical explanation of the simulation is needed
- change pdf plots to ccdf plots (emphasis on the tail of the distribution)
- Show a few regulated datasets to show that they have a high GMP-cor.

**New Panels**
- Panel A: Same as before, introduce MP distribution as correlation spectrum from Random Gaussian Matrix.
- Panel B: Still show original GMP results, emphasise that it is now the result of a random Gaussian matrix but with underlying correlations.
- Panel C: Illustrative heatmap of random sparse matrix that mimics the data structure used in the simulation, can be exaggerated to show some non-zero values in the relatively small matrix. Subplot: same heatmap but with scrambled indices within each column, add arrows to depict scrambling.
- Panel D: correlation eigenvalues of example simulation results. Use the following code to format the plot:
```python
import numpy as np
import matplotlib.pyplot as plt
import os
[pcs,pcs1] = np.load(os.path.join('..','ev_data','simulated_pcs.npy'))
# 1. Generate data
data1 = pcs[pcs>0]
data2 = pcs1[pcs1>0]
bin_width = 0.2
x1 = (1+np.sqrt(2))**2 #MP limit
x2 = np.max(pcs1)

all_data = np.concatenate([data1, data2])
bin_edges = np.arange(min(all_data), max(all_data) + bin_width, bin_width)

fig, ax = plt.subplots(figsize=(10, 6))

# 2. Plot Dataset 1 (making bars slightly thinner than the bin width)
n1, bins1, patches1 = ax.hist(data1, bins=bin_edges, width=bin_width*0.8, align='right',
                              edgecolor='black',color='#d9d9d9', alpha=0.7, density=True)

# 3. Plot Dataset 2 (shifting its bin edges forward so they sit side-by-side/overlap)
#n2, bins2, patches2 = ax.hist(data2, bins=bin_edges + (bin_width * 0.5), width=bin_width * 0.5,
#                              edgecolor='black', color='#737373', alpha=0.7, label='Dataset 2')

# 4. Color Dataset 1 (Standard solid colors)
for patch in patches1:
    bin_x = patch.get_x()
    if bin_x < x1:
        patch.set_facecolor('darkgray')
    elif x1 <= bin_x < x2:
        patch.set_facecolor('salmon')
    else:
        patch.set_facecolor('skyblue')


# 6. Add thresholds and layout elements
ax.axvline(x1, color='k', linestyle='--', alpha=0.6)
ax.axvline(x2, color='k', linestyle='--', alpha=0.6)

ax.set_xlabel(r"Eigenvalue - $\lambda$", fontsize=14)
ax.set_ylabel(r"Density - $\rho(\lambda)$", fontsize=14)
ax.set_yticks([0,0.1,0.2,0.3,0.4])
ax.grid(False)
ax.set_xlim([bin_edges[0], bin_edges[-1]+bin_width])
plt.show()
```
- Panel D subplot: histogram of `pcs1` (scrambled) directly under the above plot to compare the original to scrambled, make sure the x-axis matches so the $\lambda^{max}_{scr}$ threshold is clear in both subplots.
- Panel E: example pdf and inset ccdf of regulated dataset. use the following format:
```python
import numpy as np
import matplotlib.pyplot as plt
import os
arr = np.load(os.path.join('..','ev_data','sample_13b_filtered.npy'))
# 1. Generate data
data1 = arr[0,:]
data2 = arr[1,:]
data1 = data1[data1>0]
data2 = data2[data2>0]
alpha = 2
x1 = (1+np.sqrt(alpha))**2 #MP limit
x2 = np.max(data2)
bin_width = 0.2
all_data = np.concatenate([data1, data2])
bin_edges = np.arange(min(all_data), 10 + bin_width, bin_width)

fig, ax = plt.subplots(figsize=(10, 6))

# 2. Plot Dataset 1 (making bars slightly thinner than the bin width)
n1, bins1, patches1 = ax.hist(data1, bins=bin_edges, width=bin_width*0.5, align='left',
                              edgecolor='black',color='#d9d9d9', alpha=0.7, density=True)

# 3. Plot Dataset 2 (shifting its bin edges forward so they sit side-by-side/overlap)
n2, bins2, patches2 = ax.hist(data2, bins=bin_edges + (bin_width * 0.5), width=bin_width * 0.5,
                              edgecolor='black', color='black', alpha=0.7, label='Dataset 2',align='right', density=True)

# 4. Color Dataset 1 (Standard solid colors)
for patch in patches1:
    bin_x = patch.get_x()
    if bin_x < x2:
        patch.set_facecolor('darkgray')
    else:
        patch.set_facecolor('skyblue')


# 6. Add thresholds and layout elements
#ax.axvline(x1, color='k', linestyle='--', alpha=0.6)
#ax.axvline(x2, color='k', linestyle='--', alpha=0.6)

ax.set_xlabel(r"Eigenvalue - $\lambda$", fontsize=14)
ax.set_ylabel(r"Density - $\rho(\lambda)$", fontsize=14)
ax.set_yticks([0,0.1,0.2,0.3,0.4])
ax.grid(False)
ax.set_xlim([bin_edges[0], bin_edges[-1]+bin_width])
inset_ax = ax.inset_axes([0.4, 0.4, 0.55, 0.55])
data1 = np.sort(data1)
data2 = np.sort(data2)
p = len(data1)
cdf = np.arange(1, p + 1) / p
# 3. Calculate the CCDF (1 - CDF)
ccdf = 1 - cdf + (1/p)
noise_ind = data1<x2
inset_ax.loglog(data1[noise_ind], ccdf[noise_ind], marker='.', linestyle='-',color='darkgray', alpha=0.7, label = 'noise')
inset_ax.loglog(data1[np.invert(noise_ind)], ccdf[np.invert(noise_ind)], marker='.', linestyle='-',color='skyblue', label = 'signal')
inset_ax.loglog(data2, ccdf, marker='.', linestyle='-',color='black', alpha=0.5, label = 'scrambled')
inset_ax.set_xlim([0.1,np.max(data1)])
inset_ax.legend()
inset_ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
inset_ax.set_xlabel(r"$\lambda$", fontsize=12)
inset_ax.set_ylabel(r"CCDF", fontsize=12)
plt.show()
```
-Panel F-G: two more example datasets (use sample_2b_filtered and sample_15b_filtered in `../ev_data/`) of regulated conditions but ccdf only, use the same format as in Panel E inset.

Tile the panels to optimize the space in a full page figure, but keep order and flow of information to make sense.