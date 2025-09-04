# pynbody-config
EDGE-INFERNO configuration file and particle properties

When loading ramses output from the inferno model using pynbody, derived quantities can be imported from the inferno module for convenience.

Example:
```python
import pynbody
from inferno.derived_quantity import *
 
import matplotlib.pyplot as plt
import numpy as np

ds = pynbody.load("output_00104/")
ds.physical_units()

fig = plt.figure()
plt.xlabel("Stellar age [Gyr]")
plt.ylabel("[Fe/H]")
plt.ylim(-4.2, -0.8)
plt.xlim(11.5, 13.8)

plt.scatter(ds.s['age'].in_units('Gyr'), ds.s['FeH'], c='k', s=1, rasterized=True)

plt.show()
```
