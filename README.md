dielectric_spectroscopy_gui
=======

MATLAB scripts/functions/GUIs for reading, homogenizing and analyzing dielectric spectroscopy data collected by the Solartron platform. Optimized for analyzing ice-rich materials, but potentially adjustable to other materials.

If the format of the data input into inv_cole_gui.m (determined by load_solartron_gui.m and merge_solartron_gui.m), then inv_cole_gui.m may be used to invert dielectric spectroscopy for multiple (up to 4) Cole-Cole relaxations and DC conductivity simultaneously, with confidence bounds for each model parameter.

dielectric_spectroscopy_man_gui.pdf provides a manual for GUI operation.

The sub-directory data_example provides example data from measurements of a sample from the Vostok 5G ice core. These data can be loaded in the appropriate order to see the operation of each GUI. The order of the directories is indicative of their use in GUIs. For load_solartron_gui, the files in data_example/0_cat must be moved to a sub-directory called “cal” in your working directory prior to starting load_solartron_gui. 