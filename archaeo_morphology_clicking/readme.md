# Readme for Manzuk et al. submission to Nature

Codes are organized several scripts that call on a set of functions to perform measurements and analyses. Workflow is as follows (more details with comments of scripts and functions):

1. Use `archaeo_clicking_wrapper.m` to gather outline points of all archaeocyathid specimens from downsampled GIRI image stacks. Final manually-clicked datasets are included in the data repository under /3d_data/archaeo_clicking_matlab/.
2. Use the first two cells of `coral_ct_wrapper.m` to automatically get the branch outlines for all coral specimens and clean up spurious edges. Cleaned outline data are inclueded in the data repository under /3d_data/coral_ct_matlab/.
3. To make measurements of both archaeocyathid and coral specimens use `archaeo_measurement_wrapper.m` and `coral_ct_wrapper.m` (starting from line 60) to take the saved .mat outline data and make all morphological measurements (details in script and function comments, as well as function headers).
4. Use `stromatolite_wrapper.m` and the Cantine, et al. data (repository \misc\PreC_out_20191012.xlsx) to model stromatolites and make morphological measurements.
5. To output figures 3 and 4 stemming from morphological analysis, follow br_morphology_figures.m.

Other scripts inclueded are:
1. `coral_simulation_wrapper.m`, which takes binary tiff stacks with the synthetic growth models (repository 3d_data/coral_simulation_stacks/) and automatically traces outlines and makes morphological measurements.
2. `bahamas_reef_estimator.m`, which uses the figure from geyman et al. 2021 to estimate the map coverage of reefs on the GBB.
3. `diversity_curve_wrapper.m`, which takes the paleobiology database data (repository misc/pbdb_data.csv) to produce the lower cambrian diversity curves in figure 4.
