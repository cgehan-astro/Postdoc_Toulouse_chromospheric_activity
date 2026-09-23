### Project: Results.py

I completed this project during my postdoc at at Observatoire Midi-Pyrénées, Toulouse (2024 - 2026).


### Project overview

This project aims at providing a physical interpretation of the measurements of some physical parameters that has been measured from spectroscopic data by the pipeline presented in Activity_indices_measurements.py, namely physical indicators of the level of magnetic activity of stars that is caused by internal magnetic fields.


### Dataset

The dataset is composed of the physical parameters contained in the file Activity_indices_measurements.csv, as well as additional physical parameters from the files table_physical_parameters.dat and properties_binary_stars.txt.


### Results

Two plots are saved, representing two examples of what we can learn from the data: Evolution_oscillations_amplitude_with_MgI_activity_colored.pdf and CaHK_activity_versus_rotation_single_stars_versus_binaries.pdf.

1. The plot Evolution_oscillations_amplitude_with_MgI_activity_colored.pdf represents two parameters measured with asteroseismology: the amplitude of the excess power in the power spectrum due to stellar oscillations, as a function of the frequency of maximum oscillation power. The color code represents the level of magnetic activity of stars as probed by one the neutral magnesium spectral lines.

2. The plot CaHK_activity_versus_rotation_single_stars_versus_binaries.pdf represents the level of magnetic activity of stars as probed by the two ionised calcium H & K spectral lines, as a function of the rotation period at which stars rotate on themselves. The blue dots represent single stars, the black dots represent stars belonging to a close binary system (i.e. where the two companion stars orbit around each other on an orbital period shorter than 150 days), the red dots represent stars belonging to a wide binary system (i.e. where the two companion stars orbit around each other on an orbital period of at least 150 days).


### Key findings

 1. The plot Evolution_oscillations_amplitude_with_MgI_activity_colored.pdf indicates that, at a given frequency of maximum oscillation power that is a proxy of the degree of evolution of stars, the amplitude of stellar oscillations is anti-correlated with the level of magnetic activity, based on about 3100 measurements. This is expected because magnetic fields tend to inhibit convection in the outer envelope of stars, which triggers oscillations; a larger level of magnetic activity is therefore expected to result in a partial and even sometimes in a total suppression of oscillations.

2. The plot CaHK_activity_versus_rotation_single_stars_versus_binaries.pdf emphasizes a key result for theoretical work and simulations of tidal magnetic dynamos: stars in close binaries exhibit enhanced levels of magnetic activity compared to single stars and to stars in wide binaries. This result is unexpected because it is not explained by the fact that stars in close binaries are in fast rotation: indeed, fast-rotating single stars exhibit lower levels of magnetic acivity compared to close binaries, for similar rotation periods. Therefore, magnetic fields are amplified in stars belonging to close binaries, suggesting that tides somehow amplify magnetic fields, and the mechanisms at work are still unknown.


### Installation: with anaconda

git clone https://github.com/cgehan-astro/Postdoc_Toulouse_chromospheric_activity.git
cd Postdoc_Toulouse_chromospheric_activity
conda env create -f environment.yml
