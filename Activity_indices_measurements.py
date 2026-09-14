### Create the appropriate environment with anaconda: conda env create -f environment.yml
### python=3.9.25, ipython=8.15


### This program aims at measuring several physical parameters from spectroscopic data

import numpy as np
import pylab as plt
from astropy.io import fits
import glob
import csv
import os
import scipy.signal


### Function saving the parameters enabling the identification of the object

def coordinates_star():
	list_KIC = []
	list_ra = []
	list_dec = []
	with open('table_identify_objects.txt', 'r') as fin:
		reader=csv.reader(fin)
		stars = [[file_ind for file_ind in row] for number,row in enumerate(reader)]
		for i in range(1,np.size(stars, axis=0)):
			KIC_i = int(stars[i][0])
			ra_i = float(stars[i][2])
			dec_i = float(stars[i][3])
			list_KIC.append(KIC_i)
			list_ra.append(ra_i)
			list_dec.append(dec_i)
	all_KIC = np.array(list_KIC)
	ra_KIC = np.array(list_ra)
	dec_KIC = np.array(list_dec)
	return all_KIC, ra_KIC, dec_KIC



### Function identifying to which object corresponds each spectrum: several spectra with different filenames can correspond to the same object

def KIC_identification(dec, dec_KIC, ra, ra_KIC, all_KIC):
	arg_dist = np.sin(dec) * np.sin(dec_KIC) + np.cos(dec) * np.cos(dec_KIC) * np.cos(ra - ra_KIC)		# computing the distance between the coordinates of the object for which we have a spectrum
														# and the reference objects listed beforehand
	max_arg_dist_crit = np.where(arg_dist > 1)[0]		# dealing with values that are singular: we take the arccosine of those values later on
	min_arg_dist_crit = np.where(arg_dist < -1)[0]
	arg_dist[max_arg_dist_crit] = 1.
	arg_dist[min_arg_dist_crit] = -1.
	dist = np.arccos(arg_dist)

	dist_crit = np.where(abs(dist) == min(abs(dist)))[0]	# finding the reference object minimizing the distance between coordinates

	KIC = all_KIC[dist_crit][0].astype(int)			# associating the object for which we have a spectrum with the correct reference object
	return KIC



### Function plotting the spectrum

def plot_whole_spectrum(wavelength, flux, KIC, star_number):
	plt.figure()
	plt.plot(wavelength, flux)
	plt.xscale('log')		# plot on a logarithmic scale
	plt.xlim(min(wavelength), max(wavelength))
	plt.xlabel(r'$\lambda$' + ' (' + r'$\AA$' + ')', fontsize='x-large')
	plt.ylabel('Flux (e-)', fontsize='x-large')
	plt.title('KIC ' + str(KIC), fontsize='x-large')
	plt.savefig('whole_spectrum_KIC_' + str(KIC) + '_' + star_number + '.pdf', format='pdf')
	return


### Function selecting the wavelength range around a given spectral line suitable to measure the associated physical parameter of interest

def selecting_lambda_Doppler(lambda_ref_0, wavelength_range, wavelength, flux, inverse_variance):
	wavelength_ref = wavelength[np.where(wavelength >= wavelength_range[0])[0]]
	flux = flux[np.where(wavelength >= wavelength_range[0])[0]]
	inverse_variance = inverse_variance[np.where(wavelength >= wavelength_range[1])[0]]
	wavelength_ref = wavelength_ref[np.where(wavelength_ref <= wavelength_range[1])[0]]
	flux = flux[np.where(wavelength_ref <= wavelength_range[1])[0]]
	inverse_variance = inverse_variance[np.where(wavelength_ref <= wavelength_range[1])[0]]

	min_crit = np.where(flux == min(flux))[0]	# selecting the frequency corresponding to the local flux minimum
	lambda_ref = wavelength_ref[min_crit]
	shift = lambda_ref - lambda_ref_0		# computing the shift compared to the reference frequency
	shift = shift[0]
	return wavelength_ref, flux, inverse_variance, shift



### Function computing the percentage of unphysical negative flux around a given spectral line of interest

def compute_percentage_negative_flux(flux):
	crit_negative_flux = np.where(flux < 0)[0]
	percentage_negative_flux = float(np.size(crit_negative_flux)) / np.size(flux) * 100
	return percentage_negative_flux



### Function computing the amount of Doppler shift using only physically valid portions of the spectrum: percentage of negative flux less than 1% and SNR above 15

def compute_radial_velocity(SNR_limit, percentage_negative_flux_limit, c, percentage_negative_flux_CaK, SNR_CaK, percentage_negative_flux_CaH, SNR_CaH, percentage_negative_flux_Halpha, SNR_Halpha, percentage_negative_flux_FeI, SNR_FeI, percentage_negative_flux_CaIR, SNR_CaIR):
	if percentage_negative_flux_CaK < percentage_negative_flux_limit and np.mean(SNR_CaK) > SNR_limit:
		radial_velocity_CaK = c * shift_CaK * 10.**(-3) / lambda_ref_0_CaK
	else:
		radial_velocity_CaK =  -9999.
	if percentage_negative_flux_CaH < percentage_negative_flux_limit and np.mean(SNR_CaH) > SNR_limit:
		radial_velocity_CaH =  c * shift_CaH * 10.**(-3) / lambda_ref_0_CaH
	else:
		radial_velocity_CaH =  -9999.
	if percentage_negative_flux_Halpha < percentage_negative_flux_limit and np.mean(SNR_Halpha) > SNR_limit:
		radial_velocity_Halpha =  c * shift_Halpha * 10.**(-3) / lambda_ref_0_Halpha
	else:
		radial_velocity_Halpha =  -9999.
	if percentage_negative_flux_FeI < percentage_negative_flux_limit and np.mean(SNR_FeI) > SNR_limit:
		radial_velocity_FeI =  c * shift_FeI * 10.**(-3) / lambda_ref_0_FeI
	else:
		radial_velocity_FeI =  -9999.
	if percentage_negative_flux_CaIR < percentage_negative_flux_limit and np.mean(SNR_CaIR) > SNR_limit:
		radial_velocity_CaIR =  c * shift_CaIR * 10.**(-3) / lambda_ref_0_CaIR
	else:
		radial_velocity_CaIR =  -9999.

	radial_velocity_array = np.array([radial_velocity_CaK, radial_velocity_CaH, radial_velocity_Halpha, radial_velocity_FeI, radial_velocity_CaIR])
	radial_velocity_array = radial_velocity_array[np.where(radial_velocity_array != -9999.)[0]]
	radial_velocity = np.median(radial_velocity_array)
	return radial_velocity



### Function computing the final Doppler shift in frequency for a given spectral line of interest

def compute_lambda_shift(radial_velocity, lambda_ref_0, c, wavelength_ref, flux, dnu_triangle, width_for_lambda_shift):
	shift_lambda_ref = radial_velocity * 10.**3 * lambda_ref_0 / c		# frequency shift based on the Doppler effect
	lambda_ref = lambda_ref_0 + shift_lambda_ref				# definition of the actual frequency that is shifted compared to the reference frequency
	
	flux_equal = (np.interp(lambda_ref - dnu_triangle, wavelength_ref, flux) + np.interp(lambda_ref + dnu_triangle, wavelength_ref, flux)) / 2  # mean of the flux interpolated at the edges of the center of spectral line: target value to refine the computation of the Doppler shift
	
	## Defining a larger wavelength range around the line center
	
	lambda_ref_equal = wavelength_ref[np.where(wavelength_ref >= lambda_ref - width_for_lambda_shift * dnu_triangle)[0]]
	flux_ref_equal = flux[np.where(wavelength_ref >= lambda_ref - width_for_lambda_shift * dnu_triangle)[0]]
	lambda_ref_equal = lambda_ref_equal[np.where(lambda_ref_equal <= lambda_ref + width_for_lambda_shift * dnu_triangle)[0]]
	flux_ref_equal = flux_ref_equal[np.where(lambda_ref_equal <= lambda_ref + width_for_lambda_shift * dnu_triangle)[0]]
	
	## Separating the wavelength range around the line center: left/blue part and right/red part
	
	flux_ref_equal_1 = flux_ref_equal[np.where(lambda_ref_equal <= lambda_ref)[0]]
	lambda_ref_equal_1 = lambda_ref_equal[np.where(lambda_ref_equal <= lambda_ref)[0]]
	flux_ref_equal_2 = flux_ref_equal[np.where(lambda_ref_equal >= lambda_ref)[0]]
	lambda_ref_equal_2 = lambda_ref_equal[np.where(lambda_ref_equal >= lambda_ref)[0]]
	
	## Computation of the final refined Doppler shift in frequency: mean between the interpolated frequencies associated to the target of the flux (flux_equal defined above)
	
	lambda_ref = (np.interp(flux_equal, flux_ref_equal_1[np.argsort(flux_ref_equal_1)], lambda_ref_equal_1[np.argsort(flux_ref_equal_1)]) + np.interp(flux_equal, flux_ref_equal_2[np.argsort(flux_ref_equal_2)], lambda_ref_equal_2[np.argsort(flux_ref_equal_2)])) / 2
	return lambda_ref



### Function rebinning the spectrum using a higher resolution

def interpolated_spectrum(lambda_ref, dnu_triangle, new_R, lambda_square_C, dnu_square, wavelength, flux):
	N = int(np.round((lambda_ref + dnu_triangle - (lambda_ref - dnu_triangle)) / new_R))
	N_continuum = int(np.round((lambda_square_C + dnu_square - (lambda_square_C - dnu_square)) / new_R))
	nu_around_ref = np.linspace(lambda_ref - dnu_triangle, lambda_ref + dnu_triangle, N)	# flux surrounding the line
	flux_around_ref = np.interp(nu_around_ref, wavelength, flux)	# interpolation
	triangle_filter = scipy.signal.windows.triang(N)		# triangular filter
	nu_triangle_ref = np.linspace(lambda_ref - dnu_triangle, lambda_ref + dnu_triangle, N)
	flux_ref = np.interp(nu_triangle_ref, nu_around_ref, flux_around_ref)	# flux in the line weighted by the triangular filter
	nu_around_C = np.linspace(lambda_square_C - dnu_square, lambda_square_C + dnu_square, N_continuum)	# flux surrounding the continuum
	flux_around_C = np.interp(nu_around_C, wavelength, flux)	# interpolation
	square_filter = np.zeros((N_continuum))
	square_filter[:] = 1.		# squared filter
	nu_square_C = np.linspace(lambda_square_C - dnu_square, lambda_square_C + dnu_square, N_continuum)
	flux_C = np.interp(nu_square_C, nu_around_C, flux_around_C)		# interpolated flux in the continuum
	return triangle_filter, nu_triangle_ref, flux_ref, square_filter, nu_square_C, flux_C



### Function computating the uncertainties on the measured physical parameters
	
def compute_uncertainties(wavelength, flux, inverse_variance, lambda_ref, lambda_square_C, nu_square_C):
	crit = np.where(wavelength >= lambda_ref - (lambda_square_C - lambda_ref))[0]	# wavelength selection to compute the SNR
	wavelength = wavelength[crit]
	flux = flux[crit]
	inverse_variance = inverse_variance[crit]
	crit = np.where(wavelength <= max(nu_square_C))[0]		# wavelength selection to compute the SNR
	wavelength = wavelength[crit]
	flux = flux[crit]
	inverse_variance = inverse_variance[crit]
	SNR = flux * inverse_variance**0.5		# SNR computation
	log_sigma_S = - np.log10(np.mean(SNR)) - 0.5	# uncertainty: scaling relation with the SNR		
	sigma_S = 10.**log_sigma_S
	return wavelength, flux, inverse_variance, SNR, sigma_S



### Function plotting the spectrum centered around the H&K lines of CaII

def plot_CaHK_spectrum(wavelength, flux, KIC, star_number):
	plt.figure()
	plt.plot(wavelength_CaHK, flux_CaHK)
	plt.axvline(lambda_ref_CaK, linestyle='--', linewidth=0.5, color='k', label='Line centers')
	plt.axvline(lambda_ref_CaH, linestyle='--', linewidth=0.5, color='k')
	plt.plot(nu_triangle_ref_CaK, flux_ref_CaK, c='r', label ='Flux selection in the line centers')
	plt.plot(nu_triangle_ref_CaH, flux_ref_CaH, c='r')
	plt.plot(nu_square_B_CaHK, flux_B_CaHK, c='orange', label ='Flux selection in the nearby continuum')
	plt.plot(nu_square_R_CaHK, flux_R_CaHK, c='orange')
	plt.legend()
	plt.xlim(min(wavelength_CaHK), max(wavelength_CaHK))
	plt.xlabel(r'$\lambda$' + ' (' + r'$\AA$' + ')', fontsize='x-large')
	plt.ylabel('Flux (e-)', fontsize='x-large')
	plt.title('KIC ' + str(KIC), fontsize='x-large')
	plt.savefig('CaHK_spectrum_KIC_' + str(KIC) + '_' + star_number + '.pdf', format='pdf')
	plt.close()



### Function writing the results in an output file

def write_results(fichier, KIC, S_CaHK, sigma_S_CaHK, S_MgI, sigma_S_MgI, S_FeI, sigma_S_FeI, S_Halpha, sigma_S_Halpha, S_CaIR, sigma_S_CaIR):
	fichier.write(str(KIC) + ' ' + str('%.3f' % S_CaHK) + ' ' + str('%.3f' % sigma_S_CaHK) + ' ' + str('%.3f' % S_MgI) + ' ' + str('%.3f' % sigma_S_MgI) + ' ' + str('%.3f' % S_FeI) + ' ' + str('%.3f' % sigma_S_FeI) + ' ' + str('%.3f' % S_Halpha) + ' ' + str('%.3f' % sigma_S_Halpha) + ' ' + str('%.3f' % S_CaIR) + ' ' + str('%.3f' % sigma_S_CaIR) + '\n')
	return
	
	

## List of spectra to analyse 

path = './spec*.fits.gz'
files = glob.glob(path)
files_number = np.size(files, axis=0)


## Defining the output file to save the results

output_path = "./output_measurements.txt"
if os.path.exists(output_path) == False:
	fichier = open(output_path, "w")
if os.path.exists(output_path) == True:
	fichier = open(output_path, "a")	# append the results if the file already exists to not erase previous results


## Saving the parameters enabling the identification of the object

all_KIC, ra_KIC, dec_KIC = coordinates_star()


## Definition of the frequencies of interest  around each spectral line for each physical parameter to measure

lambda_ref_0_Halpha = 6562.801
Halpha_range = np.array([6530, 6590])

lambda_ref_0_CaK = 3933.664
CaK_range = np.array([3850, 3950])

lambda_ref_0_CaH = 3968.470
CaH_range = np.array([3951, 4025])

lambda_ref_0_FeI = 6495
FeI_range = np.array([6475, 6524])

lambda_ref_0_CaIR = 8542
CaIR_range = np.array([8510, 8600])

lambda_ref_0_MgI = 5184
MgI_range = np.array([5130, 5250])

lambda_square_B_0_CaHK = 3901.070
lambda_square_R_0_CaHK = 4001.070
lambda_square_C_0_Halpha = 6605.
lambda_square_C_0_FeI = 6605.
lambda_square_C_0_MgI = 5240.
lambda_square_C_0_CaIR = 8620.

dnu_triangle_CaHK = 2.18				# FWHM of the triangle function
dnu_triangle_Halpha = 1.5				# FWHM of the triangle function
dnu_triangle_FeI = 1.5				# FWHM of the triangle function
dnu_triangle_MgI = 1.5				# FWHM of the triangle function
dnu_triangle_CaIR = 1.5				# FWHM of the triangle function

dnu_square_CaHK = 10
dnu_square_Halpha = 10
dnu_square_FeI = 10
dnu_square_MgI = 10
dnu_square_CaIR = 10

width_for_lambda_shift = 5.


## Definition of some constants

SNR_limit = 15.		# SNR threshold above which the signal is significant compared to the noise
SNR_limit_final = 8.	# final SNR threshold above which the signal is significant compared to the noise
percentage_negative_flux_limit = 1.	# limit on the percentage of unphysical negative flux above which the measurements are not valid
c = 3.*10.**8		# light speed
new_R = 10.**(-5)	# new, higher, spectral resolution


## Beginning of the main program

for i in range(files_number):		# loop on the different spectra

	hdul = fits.open(files[i])
	
	ra = hdul[0].header['RA']
	dec = hdul[0].header['DEC']

	KIC = KIC_identification(dec, dec_KIC, ra, ra_KIC, all_KIC) 		# identifying to which object corresponds each spectrum
	
	
	## Reading the physical parameters needed to exploit the spectrum
	
	data = hdul[0].data
	wavelength = data[2,:]
	flux = data[0,:]
	inverse_variance = data[1,:]	# used to estimate the noise level
	
	
	## Plot the spectrum
	
	star_number = files[i].replace('./spec-', '')
	star_number = star_number.replace('.fits.gz', '')
	plot_whole_spectrum(wavelength, flux, KIC, star_number)
	

	## Selecting the wavelength ranges around each spectral line suitable to measure each desired physical parameter

	wavelength_ref_Halpha, flux_Halpha, inverse_variance_H_alpha, shift_Halpha = selecting_lambda_Doppler(lambda_ref_0_Halpha, Halpha_range, wavelength, flux, inverse_variance)
	wavelength_ref_CaK, flux_CaK, inverse_variance_CaK, shift_CaK = selecting_lambda_Doppler(lambda_ref_0_CaK, CaK_range, wavelength, flux, inverse_variance)
	wavelength_ref_CaH, flux_CaH, inverse_variance_CaH, shift_CaH = selecting_lambda_Doppler(lambda_ref_0_CaH, CaH_range, wavelength, flux, inverse_variance)
	wavelength_ref_FeI, flux_FeI, inverse_variance_FeI, shift_FeI = selecting_lambda_Doppler(lambda_ref_0_FeI, FeI_range, wavelength, flux, inverse_variance)
	wavelength_ref_CaIR, flux_CaIR, inverse_variance_CaIR, shift_CaIR = selecting_lambda_Doppler(lambda_ref_0_CaIR, CaIR_range, wavelength, flux, inverse_variance)
	wavelength_ref_MgI, flux_MgI, inverse_variance_MgI, shift_MgI = selecting_lambda_Doppler(lambda_ref_0_MgI, MgI_range, wavelength, flux, inverse_variance)


	## Computing the percentage of unphysical negative flux around each spectral line of interest

	percentage_negative_flux_Halpha = compute_percentage_negative_flux(flux_Halpha)
	percentage_negative_flux_CaK = compute_percentage_negative_flux(flux_CaK)
	percentage_negative_flux_CaH = compute_percentage_negative_flux(flux_CaH)
	percentage_negative_flux_FeI = compute_percentage_negative_flux(flux_FeI)
	percentage_negative_flux_CaIR = compute_percentage_negative_flux(flux_CaIR)


	## Estimating the SNR around each spectral line of interest

	SNR_Halpha = flux_Halpha * inverse_variance_H_alpha**0.5
	SNR_CaK = flux_CaK * inverse_variance_CaK**0.5
	SNR_CaH = flux_CaH * inverse_variance_CaH**0.5
	SNR_FeI = flux_FeI * inverse_variance_FeI**0.5
	SNR_CaIR = flux_CaIR * inverse_variance_CaIR**0.5


	## Measuring the amount of Doppler shift using only physically valid portions of the spectrum: percentage of negative flux less than 1% and SNR above 15

	radial_velocity = compute_radial_velocity(SNR_limit, percentage_negative_flux_limit, c, percentage_negative_flux_CaK, SNR_CaK, percentage_negative_flux_CaH, SNR_CaH, percentage_negative_flux_Halpha, SNR_Halpha, percentage_negative_flux_FeI, SNR_FeI, percentage_negative_flux_CaIR, SNR_CaIR)


	## Computing the Doppler shift in frequency for the spectral lines of interest and the continnuum nearby
	
	lambda_ref_Halpha = compute_lambda_shift(radial_velocity, lambda_ref_0_Halpha, c, wavelength_ref_Halpha, flux_Halpha, dnu_triangle_Halpha, width_for_lambda_shift)
	shift_lambda_C_Halpha = radial_velocity * 10.**3 * lambda_square_C_0_Halpha / c
	lambda_square_C_Halpha = lambda_square_C_0_Halpha + shift_lambda_C_Halpha

	lambda_ref_FeI = compute_lambda_shift(radial_velocity, lambda_ref_0_FeI, c, wavelength_ref_FeI, flux_FeI, dnu_triangle_FeI, width_for_lambda_shift)
	shift_lambda_C_FeI = radial_velocity * 10.**3 * lambda_square_C_0_FeI / c
	lambda_square_C_FeI = lambda_square_C_0_FeI + shift_lambda_C_FeI

	lambda_ref_CaIR = compute_lambda_shift(radial_velocity, lambda_ref_0_CaIR, c, wavelength_ref_CaIR, flux_CaIR, dnu_triangle_CaIR, width_for_lambda_shift)
	shift_lambda_C_CaIR = radial_velocity * 10.**3 * lambda_square_C_0_CaIR / c
	lambda_square_C_CaIR = lambda_square_C_0_CaIR + shift_lambda_C_CaIR

	lambda_ref_MgI = compute_lambda_shift(radial_velocity, lambda_ref_0_MgI, c, wavelength_ref_MgI, flux_MgI, dnu_triangle_MgI, width_for_lambda_shift)
	shift_lambda_C_MgI = radial_velocity * 10.**3 * lambda_square_C_0_MgI / c
	lambda_square_C_MgI = lambda_square_C_0_MgI + shift_lambda_C_MgI
	
	lambda_ref_CaK = compute_lambda_shift(radial_velocity, lambda_ref_0_CaK, c, wavelength_ref_CaK, flux_CaK, dnu_triangle_CaHK, width_for_lambda_shift)
	lambda_ref_CaH = compute_lambda_shift(radial_velocity, lambda_ref_0_CaH, c, wavelength_ref_CaH, flux_CaH, dnu_triangle_CaHK, width_for_lambda_shift)
	shift_lambda_B_CaHK = radial_velocity * 10. ** 3 * lambda_square_B_0_CaHK / c
	lambda_square_B_CaHK = lambda_square_B_0_CaHK + shift_lambda_B_CaHK
	shift_lambda_R_CaHK = radial_velocity * 10. ** 3 * lambda_square_R_0_CaHK / c
	lambda_square_R_CaHK = lambda_square_R_0_CaHK + shift_lambda_R_CaHK


	## Rebinning the spectrum using a higher resolution: more accurate integration later on by interpolating to not underestimate the flux in a given bandpass
	
	R = new_R / np.mean(np.diff(wavelength))	# Ratio of the initial resolution and the new one
	
	triangle_filter_Halpha, nu_triangle_ref_Halpha, flux_ref_Halpha, square_filter_Halpha, nu_square_C_Halpha, flux_C_Halpha = interpolated_spectrum(lambda_ref_Halpha, dnu_triangle_Halpha, new_R, lambda_square_C_Halpha, dnu_square_Halpha, wavelength, flux)
	
	triangle_filter_FeI, nu_triangle_ref_FeI, flux_ref_FeI, square_filter_FeI, nu_square_C_FeI, flux_C_FeI = interpolated_spectrum(lambda_ref_FeI, dnu_triangle_FeI, new_R, lambda_square_C_FeI, dnu_square_FeI, wavelength, flux)
	
	triangle_filter_CaIR, nu_triangle_ref_CaIR, flux_ref_CaIR, square_filter_CaIR, nu_square_C_CaIR, flux_C_CaIR = interpolated_spectrum(lambda_ref_CaIR, dnu_triangle_CaIR, new_R, lambda_square_C_CaIR, dnu_square_CaIR, wavelength, flux)
	
	triangle_filter_MgI, nu_triangle_ref_MgI, flux_ref_MgI, square_filter_MgI, nu_square_C_MgI, flux_C_MgI = interpolated_spectrum(lambda_ref_MgI, dnu_triangle_MgI, new_R, lambda_square_C_MgI, dnu_square_MgI, wavelength, flux)
	
	N_CaHK = int(np.round((lambda_ref_CaK + dnu_triangle_CaHK - (lambda_ref_CaK - dnu_triangle_CaHK)) / new_R))
	N_continuum_CaHK = int(np.round((lambda_square_B_CaHK + dnu_square_CaHK - (lambda_square_B_CaHK - dnu_square_CaHK)) / new_R))
	nu_around_ref_CaK = np.linspace(lambda_ref_CaK - dnu_triangle_CaHK, lambda_ref_CaK + dnu_triangle_CaHK, N_CaHK)		# flux surrounding the line
	nu_around_ref_CaH = np.linspace(lambda_ref_CaH - dnu_triangle_CaHK, lambda_ref_CaH + dnu_triangle_CaHK, N_CaHK) 	# flux surrounding the line
	flux_around_ref_CaK = np.interp(nu_around_ref_CaK, wavelength, flux)		# interpolation
	flux_around_ref_CaH = np.interp(nu_around_ref_CaH, wavelength, flux)		# interpolation
	nu_around_B_CaHK = np.linspace(lambda_square_B_CaHK - dnu_square_CaHK, lambda_square_B_CaHK + dnu_square_CaHK, N_continuum_CaHK)	# flux surrounding the continuum
	nu_around_R_CaHK = np.linspace(lambda_square_R_CaHK - dnu_square_CaHK, lambda_square_R_CaHK + dnu_square_CaHK, N_continuum_CaHK)	# flux surrounding the continuum
	flux_around_B_CaHK = np.interp(nu_around_B_CaHK, wavelength, flux)		# interpolation
	flux_around_R_CaHK = np.interp(nu_around_R_CaHK, wavelength, flux)		# interpolation
	triangle_filter_CaHK = scipy.signal.windows.triang(N_CaHK)			# triangular filter
	nu_triangle_ref_CaK = np.linspace(lambda_ref_CaK - dnu_triangle_CaHK, lambda_ref_CaK + dnu_triangle_CaHK, N_CaHK)
	flux_ref_CaK = np.interp(nu_triangle_ref_CaK, nu_around_ref_CaK, flux_around_ref_CaK)		# flux in the line weighted by the triangular filter
	nu_triangle_ref_CaH = np.linspace(lambda_ref_CaH - dnu_triangle_CaHK, lambda_ref_CaH + dnu_triangle_CaHK, N_CaHK)
	flux_ref_CaH = np.interp(nu_triangle_ref_CaH, nu_around_ref_CaH, flux_around_ref_CaH)		# flux in the line weighted by the triangular filter
	square_filter_CaHK = np.zeros((N_continuum_CaHK))
	square_filter_CaHK[:] = 1.		# squared filter
	nu_square_B_CaHK = np.linspace(lambda_square_B_CaHK - dnu_square_CaHK, lambda_square_B_CaHK + dnu_square_CaHK, N_continuum_CaHK)
	flux_B_CaHK = np.interp(nu_square_B_CaHK, nu_around_B_CaHK, flux_around_B_CaHK)		# interpolated flux in the continuum
	nu_square_R_CaHK = np.linspace(lambda_square_R_CaHK - dnu_square_CaHK, lambda_square_R_CaHK + dnu_square_CaHK, N_continuum_CaHK)
	flux_R_CaHK = np.interp(nu_square_R_CaHK, nu_around_R_CaHK, flux_around_R_CaHK)		# interpolated flux in the continuum

	
	## Scaled flux in the line and continuum bands
	
	F_lambda_ref_Halpha = flux_ref_Halpha * R
	F_lambda_C_Halpha = flux_C_Halpha * R
	
	F_lambda_ref_FeI = flux_ref_FeI * R
	F_lambda_C_FeI = flux_C_FeI * R
	
	F_lambda_ref_CaIR = flux_ref_CaIR * R
	F_lambda_C_CaIR = flux_C_CaIR * R
	
	F_lambda_ref_MgI = flux_ref_MgI * R
	F_lambda_C_MgI = flux_C_MgI * R
	
	F_lambda_ref_CaK = flux_ref_CaK * R
	F_lambda_ref_CaH = flux_ref_CaH * R
	F_lambda_B_CaHK = flux_B_CaHK * R
	F_lambda_R_CaHK = flux_R_CaHK * R

	
	## Integrated flux in the line and continuum bands
	
	F_Delta_ref_Halpha = 1./(max(nu_triangle_ref_Halpha) - min(nu_triangle_ref_Halpha)) * np.sum(F_lambda_ref_Halpha * triangle_filter_Halpha)
	F_Delta_C_Halpha = 1./(max(nu_square_C_Halpha) - min(nu_square_C_Halpha)) * np.sum(F_lambda_C_Halpha * square_filter_Halpha)
	
	F_Delta_ref_FeI = 1./(max(nu_triangle_ref_FeI) - min(nu_triangle_ref_FeI)) * np.sum(F_lambda_ref_FeI * triangle_filter_FeI)
	F_Delta_C_FeI = 1./(max(nu_square_C_FeI) - min(nu_square_C_FeI)) * np.sum(F_lambda_C_FeI * square_filter_FeI)
	
	F_Delta_ref_CaIR = 1./(max(nu_triangle_ref_CaIR) - min(nu_triangle_ref_CaIR)) * np.sum(F_lambda_ref_CaIR * triangle_filter_CaIR)
	F_Delta_C_CaIR = 1./(max(nu_square_C_CaIR) - min(nu_square_C_CaIR)) * np.sum(F_lambda_C_CaIR * square_filter_CaIR)
	
	F_Delta_ref_MgI = 1./(max(nu_triangle_ref_MgI) - min(nu_triangle_ref_MgI)) * np.sum(F_lambda_ref_MgI * triangle_filter_MgI)
	F_Delta_C_MgI = 1./(max(nu_square_C_MgI) - min(nu_square_C_MgI)) * np.sum(F_lambda_C_MgI * square_filter_MgI)
	
	F_Delta_ref_CaK = 1. / (max(nu_triangle_ref_CaK) - min(nu_triangle_ref_CaK)) * np.sum(F_lambda_ref_CaK * triangle_filter_CaHK)
	F_Delta_ref_CaH = 1. / (max(nu_triangle_ref_CaH) - min(nu_triangle_ref_CaH)) * np.sum(F_lambda_ref_CaH * triangle_filter_CaHK)
	F_Delta_B_CaHK = 1. / (max(nu_square_B_CaHK) - min(nu_square_B_CaHK)) * np.sum(F_lambda_B_CaHK * square_filter_CaHK)
	F_Delta_R_CaHK = 1. / (max(nu_square_R_CaHK) - min(nu_square_R_CaHK)) * np.sum(F_lambda_R_CaHK * square_filter_CaHK)
	
	
	## Measurement of the physical parameters of interest: activity indices
	
	S_Halpha = F_Delta_ref_Halpha / F_Delta_C_Halpha
	S_FeI = F_Delta_ref_FeI / F_Delta_C_FeI
	S_CaIR = F_Delta_ref_CaIR / F_Delta_C_CaIR
	S_MgI = F_Delta_ref_MgI / F_Delta_C_MgI
	S_CaHK = (F_Delta_ref_CaH + F_Delta_ref_CaK) / (F_Delta_B_CaHK + F_Delta_R_CaHK)
	

	## Computation of the uncertainties on the measured physical parameters
	
	wavelength_Halpha, flux_Halpha, inverse_variance_Halpha, SNR_Halpha, sigma_S_Halpha = compute_uncertainties(wavelength, flux, inverse_variance, lambda_ref_Halpha, lambda_square_C_Halpha, nu_square_C_Halpha)
	
	wavelength_FeI, flux_FeI, inverse_variance_FeI, SNR_FeI, sigma_S_FeI = compute_uncertainties(wavelength, flux, inverse_variance, lambda_ref_FeI, lambda_square_C_FeI, nu_square_C_FeI)
	
	wavelength_CaIR, flux_CaIR, inverse_variance_CaIR, SNR_CaIR, sigma_S_CaIR = compute_uncertainties(wavelength, flux, inverse_variance, lambda_ref_CaIR, lambda_square_C_CaIR, nu_square_C_CaIR)
	
	wavelength_MgI, flux_MgI, inverse_variance_MgI, SNR_MgI, sigma_S_MgI = compute_uncertainties(wavelength, flux, inverse_variance, lambda_ref_MgI, lambda_square_C_MgI, nu_square_C_MgI)

	crit_CaHK = np.where(wavelength >= min(nu_square_B_CaHK))[0]			# wavelength selection to compute the SNR
	wavelength_CaHK = wavelength[crit_CaHK]
	flux_CaHK = flux[crit_CaHK]
	inverse_variance_CaHK = inverse_variance[crit_CaHK]
	crit_CaHK = np.where(wavelength_CaHK <= max(nu_square_R_CaHK))[0]		# wavelength selection to compute the SNR
	wavelength_CaHK = wavelength_CaHK[crit_CaHK]
	flux_CaHK = flux_CaHK[crit_CaHK]
	inverse_variance_CaHK = inverse_variance_CaHK[crit_CaHK]
	SNR_CaHK = flux_CaHK * inverse_variance_CaHK**0.5			# SNR computation
	log_sigma_S_CaHK = - np.log10(np.mean(SNR_CaHK)) - 0.5  		# uncertainty: scaling relation with the SNR	
	sigma_S_CaHK = 10.**log_sigma_S_CaHK


	## Final check on the physical validity of the measurements: valid only if the final percentage of negative flux is below 1% and the final SNR is above 8

	crit_negative_flux_Halpha = np.where(flux_Halpha < 0)[0]
	percentage_negative_flux_Halpha = float(np.size(crit_negative_flux_Halpha)) / np.size(flux_Halpha) * 100
	if percentage_negative_flux_Halpha >= percentage_negative_flux_limit and np.mean(SNR_Halpha) <= SNR_limit_final:
		S_Halpha = -9999.
		sigma_S_Halpha = -9999.

	crit_negative_flux_FeI = np.where(flux_FeI < 0)[0]
	percentage_negative_flux_FeI = float(np.size(crit_negative_flux_FeI)) / np.size(flux_FeI) * 100
	if percentage_negative_flux_FeI >= percentage_negative_flux_limit and np.mean(SNR_FeI) <= SNR_limit_final:
		S_FeI = -9999.
		sigma_S_FeI = -9999.

	crit_negative_flux_CaIR = np.where(flux_CaIR < 0)[0]
	percentage_negative_flux_CaIR = float(np.size(crit_negative_flux_CaIR)) / np.size(flux_CaIR) * 100
	if percentage_negative_flux_CaIR >= percentage_negative_flux_limit and np.mean(SNR_CaIR) <= SNR_limit_final:
		S_CaIR = -9999.
		sigma_S_CaIR = -9999.
		
	crit_negative_flux_MgI = np.where(flux_MgI < 0)[0]
	percentage_negative_flux_MgI = float(np.size(crit_negative_flux_MgI)) / np.size(flux_MgI) * 100
	if percentage_negative_flux_MgI >= percentage_negative_flux_limit and np.mean(SNR_MgI) <= SNR_limit_final:
		S_MgI = -9999.
		sigma_S_MgI = -9999.
		
	crit_negative_flux_CaHK = np.where(flux_CaHK < 0)[0]
	percentage_negative_flux_CaHK = float(np.size(crit_negative_flux_CaHK)) / np.size(flux_CaHK) * 100
	if percentage_negative_flux_CaHK >= percentage_negative_flux_limit and np.mean(SNR_CaHK) <= SNR_limit_final:
		S_CaHK = -9999.
		sigma_S_CaHK = -9999.


	## Plot the spectrum centered around the H&K lines of CaII
	
	plot_CaHK_spectrum(wavelength, flux, KIC, star_number)


	## Saving the results in an output file

	write_results(fichier, KIC, S_CaHK, sigma_S_CaHK, S_MgI, sigma_S_MgI, S_FeI, sigma_S_FeI, S_Halpha, sigma_S_Halpha, S_CaIR, sigma_S_CaIR)

fichier.close()
print(str(i+1) + '. KIC ' + str(KIC) + '\n')

