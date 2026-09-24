### Project: Activity_indices_measurements.py

I completed this project during my postdoc at Observatoire Midi-Pyrénées, Toulouse (2024 - 2026).


### Project overview

This project aims at measuring several physical parameters from spectroscopic data, which are physical indicators of the level of magnetic activity of stars that is caused by internal magnetic fields. Magnetic fields generate a non-thermal heating for some stellar spectral lines that are in absorption, which results in excess flux at the centre of those absorption lines. The depth of such absorption spectral lines is therefore a proxy of the strength of stellar magnetic fields, and measuring this depth provides us with a measurement of the level of magnetic activity. Understanding the manifestations of magnetic fields at the surface of stars and in their atmosphere is crucial for exoplanetary science, because magnetic activity hampers the detection of exoplanets and impacts their habitability.


### Dataset

The dataset is composed of low-resolution (resolving power of 1800) optical spectra obtained by the ground-based spectrograph LAMOST, i.e. the flux coming from stars decomposed as a function of the wavelength; here, only one spectrum is provided as an example: spec-*.fits.gz. 


### Methodology

## Data preprocessing

There are two data preprocessing steps. The first step consists in identifying to which star corresponds each spectrum, because one star can have several spectra in the database. The correct identification is based on the file table_identify_objects.txt, which contains the identifiers of stars and their celestial coordinates; the correct star associated to a given spectrum is identified by finding the identifier corresponding to the minimum distance between its coordinate and the coordinates of the star that is analyzed. The second step consists in correcting the spectrum from the radial velocity of the star, i.e. the Doppler effect due to the proper motion of the star along the line-of-sight, which shifts all the wavelengths (towards the blue if the star is getting closer from us, towards the red if the star is getting further from us). The radial velocity of the star and associated frequency shift to apply are computed based on the difference in the wavelengths of the H-alpha line observed in the spectrum and in the rest frame (in the absence of Doppler effect).


## Model

The model is composed of a set of bandpasses by which the flux is mutliplied: triangular bandpasses centered on the center of the absorption lines that are sensitive to magnetic activity, and square bandpasses in the continuum nearby.


## Measurement

Each indicator of the activity level is computed as the ratio between the sum of the integrated fluxes in the center of the absorption lines multiplied by the triangular bandpasses, and the sum of the integrated fluxes in the continuum nearby multiplied by the square bandpasses. The uncertainties are computed through a scaling relation with the signal-to-noise ratio, which depends on the flux. When several spectra are available for a given star, the consolidated indicator of the activity level corresponds to the median of the measured values, and its consolidated uncertainty corresponds to the mean square-root of the sum of the squared uncertainties. Since the correction of background light can result in a nonphysical negative flux, only the measurements associated with a high-enough signal-to-noise ratio combined with a negligible fraction of negative flux near the spectral line of interest are considered valid and kept.


### Results

The results are saved in the file output_measurements.txt; the first column corresponds to the identifier of the star analyzed while the other columns corresponds to 5 indicators of the activity level and their associated uncertainties, as probed by the ionised calcium lines H (3968.470 Å in the rest frame) & K (3933.664 Å in the rest frame), one of the neutral magnesium lines (5184 Å in the rest frame), one the neutral iron lines (6495 Å in the rest frame), the H-alpha line (6562.801 Å in the rest frame) and one the ionised infrared calcium lines (8542 Å in the rest frame). The whole optical spectrum is saved as whole_spectrum_KIC_*.pdf and shows the numerous absorption spectral lines as well; the overall spectrum follows the curved shape of a black body. The spectrum centered around the ionised calcium lines H & K is saved as CaHK_spectrum_KIC_*.pdf: the line centers are identified by the vertical black dashed lines while the flux ranges used to compute the indicator of the activity level are represented in red in the line centers and in orange in the continuum nearby; the effect of magnetic fields are cleraly visible under the form of excess emission in the center of the absorption lines.


### Key findings

This project led to the publication of a scientific article in the international review Astronomy and Astrophysics in 2018: Gehan et al. 2022, A&A, vol. 668, A116. The results presented in this paper are central to the astrophysics community, in particular for theoretical work and simulations of tidal magnetic dynamos. 


### Installation: with anaconda

git clone https://github.com/cgehan-astro/Postdoc_Toulouse_chromospheric_activity.git
cd Postdoc_Toulouse_chromospheric_activity
conda env create -f environment.yml
