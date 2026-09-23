### Organisation of the repository

There are 3 programs in this repository.


### 1. Activity_indices_measurements.py: description in readme_Activity_indices_measurements.txt

This program takes as an input the file spec-*.txt and produces as outputs the file output_measurements.txt as well as the plots CaHK_spectrum_KIC_*.pdf and whole_spectrum_KIC_*.pdf.


### 2. Results.py: description in readme_Results.txt

This program takes as an input the files Activity_indices_measurements.csv, table_physical_parameters.dat and properties_binary_stars.txt, and produces as outputs the plots Evolution_oscillations_amplitude_with_MgI_activity_colored.pdf and CaHK_activity_versus_rotation_single_stars_versus_binaries.pdf.


### 3. Pie_chart_classification.py: description in readme_Pie_chart_classification.txt

This program produces as an output the plot Pie_chart_classification.pdf.


### Installation: with anaconda

git clone https://github.com/cgehan-astro/Postdoc_Toulouse_chromospheric_activity.git
cd Postdoc_Toulouse_chromospheric_activity
conda env create -f environment.yml
