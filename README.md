# detection and analysis of mucus halo around marine snow
This analysis follows the PIV analysis given in the PIV-gravity-machine-on-oceans repository. See the dryrad dataset for further details: Chajwa, Rahul et al. (2024). Hidden comet-tails of marine snow impede ocean-based carbon sequestration [Dataset]. Dryad. https://doi.org/10.5061/dryad.v15dv4253. The code sequentially goes through each dataset and displays the PIV flow image, and prompts user for input identifier '0' or '1' for bad or good dataset respectively. This code evaluates the hydrodynamic radius (no-flow region) around marine snow through PIV data (as seen below)
![mucus_halo](https://github.com/user-attachments/assets/65e518b9-2e4e-44e7-a697-898aadbf81dd)

The red dotted ellipse is fitted around region of zero velocity using thresholding the vertical velocity profile. 

The code further yields both PIV image with visible particle and flow characteristics.

The partice shape, size and flow measure resulting from this analysis is explained in https://doi.org/10.48550/arXiv.2310.01982, which can be plotted separately. 
