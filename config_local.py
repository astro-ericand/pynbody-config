import os
base = os.environ.get("TANGOS_SIMULATION_FOLDER", "./databases/")

property_modules = os.environ.get("TANGOS_PROPERTY_MODULES", "edge_star_by_star_tangos_properties")
property_modules = property_modules.split(",")
property_modules = map(str.strip, property_modules)

# No thinning on merger trees
mergertree_min_fractional_weight = 0.00 # as a fraction of the weight of the strongest link from each halo
mergertree_min_fractional_NDM = 0.00 # as a fraction of the most massive halo at each timestep - set to zero for no thinning
mergertree_max_nhalos = 50 # maximum number of halos per step - discard the least massive ones

# Increase time to build merger trees considering no thinning
mergertree_timeout = 100.0 # seconds before abandoning the construction of a merger tree in the web interface

# Increase number of timesteps since these runs are actually very time-fine
mergertree_max_hops = 1000 # maximum number of timesteps to scan
num_multihops_max_default = 1000     # the maximum number of links to follow when searching for related halos

# Decrease relative time since these runs are very fine in time
max_relative_time_difference = 1e-5     # the maximum fractional difference in time between two contemporaneous timesteps when searching for related halos


