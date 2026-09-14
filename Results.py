### Create the appropriate environment with anaconda: conda env create -f environment.yml
### python=3.9.25, ipython=8.15


import numpy as np
import pylab as plt
import csv
import matplotlib as mpl


### Retrieving the raw the data

file = open('Activity_indices_measurements.csv')
csvreader = csv.reader(file)
rows = []
for row in csvreader:
	rows.append(row)
data_ini = np.array(rows)[1:,(0,10,11,12,13)]
list_data = []
for i in range(np.size(data_ini, axis=0)):
    if list(data_ini[i,1]):
        list_data.append(data_ini[i,:])
data = np.array(list_data).astype(float)
KIC = data[:,0]			# object identifier
SCaHK = data[:,1]		# CaHK activity indicator
sigma_SCaHK = data[:,2]		# uncertainty on the CaHK activity indicator
SMgI = data[:,3]		# MgI activity indicator
sigma_SMgI = data[:,4]		# uncertainty on the MgI activity indicator



### Crossmatching the data with another table of parameters: obtaining the rotation period, the oscillations amplitude and a proxy of evolution

data_params = np.loadtxt('table_physical_parameters.dat', usecols=(0,7,11,22))		# retrieving the parameters of interest
KIC_params = data_params[:,0]
nu_max_params = data_params[:,1]
Hmax_params = data_params[:,2]
Prot_params = data_params[:,3]
list_Prot = []
list_nu_max = []
list_Hmax = []
for i in range(KIC.size):					# associating the parameters with the correct object from the raw data
	crit = np.where(KIC_params == KIC[i])[0]		# finding the object identifier in the parameters file that corresponds to the object identifier in the raw data
	if list(crit):
		list_nu_max.append(nu_max_params[crit])
		list_Hmax.append(Hmax_params[crit])
		list_Prot.append(Prot_params[crit])
Prot = np.array(list_Prot)[:,0]		# rotation period
nu_max = np.array(list_nu_max)		# proxy of evolution
Hmax = np.array(list_Hmax)		# amplitude of oscillations



### Selecting data with valid parameters: positive nu_max, Hmax, and SMgI

crit_valid = np.where(nu_max > 0)[0]		# valid nu_max
KIC_valid = KIC[crit_valid]
nu_max_valid = nu_max[crit_valid]
Hmax_valid = Hmax[crit_valid]
SMgI_valid = SMgI[crit_valid]
sigma_SMgI_valid = sigma_SMgI[crit_valid]
crit_valid = np.where(Hmax_valid > 0)[0]	# valid Hmax
nu_max_valid = nu_max_valid[crit_valid]
Hmax_valid = Hmax_valid[crit_valid]
SMgI_valid = SMgI_valid[crit_valid]
sigma_SMgI_valid = sigma_SMgI_valid[crit_valid]
crit_valid = np.where(SMgI_valid > 0)[0]	# valid SMgI
nu_max_valid = nu_max_valid[crit_valid]
Hmax_valid = Hmax_valid[crit_valid]
SMgI_valid = SMgI_valid[crit_valid]
sigma_SMgI_valid = sigma_SMgI_valid[crit_valid]



### Plotting Hmax versus nu_max with a color code representing SMgI

cm = mpl.cm.jet
fig, ax = plt.subplots()
cax = fig.add_axes([0.55, 0.8, 0.3, 0.04])
sc = ax.scatter(nu_max_valid, Hmax_valid, c=SMgI_valid, s=20, vmin=0.45, vmax=0.82, cmap=cm, marker='o')
fig.colorbar(sc, cax=cax, ticks=[0.45, 0.54, 0.63, 0.72, 0.82], orientation='horizontal', label='$S_{MgI}$')
ax.set_xlabel(r'$\nu_{max}$' + ' (' r'$\mu$' + 'Hz)', fontsize='x-large')
ax.set_ylabel(r'$H_{max}$' + ' (ppm' + r'$^2$' + ' ' + r'$\mu$' + 'Hz' + r'$^{-1}$' + ')', fontsize='x-large')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(10,300)
ax.set_ylim(3,50000)
plt.savefig('./Evolution_oscillations_amplitude_with_MgI_activity_colored.pdf', format='pdf')
plt.close()



### Selecting data with other valid parameters: positive Prot and SCaHK

crit_valid = np.where(Prot > 0)[0]		# valid Prot
KIC_valid = KIC[crit_valid]
Prot_valid = Prot[crit_valid]
SCaHK_valid = SCaHK[crit_valid]
sigma_SCaHK_valid = sigma_SCaHK[crit_valid]
crit_valid = np.where(SCaHK_valid > 0)[0]	# valid SCaHK
KIC_valid = KIC_valid[crit_valid]
Prot_valid = Prot_valid[crit_valid]
SCaHK_valid = SCaHK_valid[crit_valid]
sigma_SCaHK_valid = sigma_SCaHK_valid[crit_valid]



### Crossmatching the data with another table of parameters: identifying stars that belong to a binary system

data_binaries = np.loadtxt('properties_binary_stars.txt', skiprows=1)		# retrieving the parameters of interest
KIC_binaries_all = data_binaries[:,0]
Prot_binaries_all = data_binaries[:,3]
Porb_binaries_all = data_binaries[:,4]
crit_binaries_valid = np.where(Porb_binaries_all > 0.)[0]	# selecting valid binaries: positive orbital period
KIC_binaries_valid = KIC_binaries_all[crit_binaries_valid]
Prot_binaries_valid = Prot_binaries_all[crit_binaries_valid]
Porb_binaries_valid = Porb_binaries_all[crit_binaries_valid]
list_crit_binaries = []
list_Porb_binaries = []
for i in range(KIC_binaries_valid.size):					# associating the parameters with the correct object from the raw data
	crit = np.where(KIC_valid == KIC_binaries_valid[i])[0]		# finding the object identifier in the binaries file that corresponds to the object identifier in the data
	if list(crit):
		list_crit_binaries.append(crit[0])
		list_Porb_binaries.append(Porb_binaries_valid[i])
crit_binaries = np.array(list_crit_binaries)
KIC_binaries = KIC_valid[crit_binaries]
Porb_binaries = np.array(list_Porb_binaries)
Prot_binaries = Prot_valid[crit_binaries]
SCaHK_binaries = SCaHK_valid[crit_binaries]
sigma_SCaHK_binaries = sigma_SCaHK_valid[crit_binaries]
crit_close_binaries = np.where(Porb_binaries <= 150)[0]		# discriminating close (orbital period equal to or below 150 days) and wide (orbital period above 150 days) binaries
crit_wide_binaries = np.where(Porb_binaries > 150)[0]



### Plotting SCaHK versus Prot by distinguishing single stars, close binaries and wide binaries

plt.figure()
plt.errorbar(Prot_valid, SCaHK_valid, yerr=sigma_SCaHK_valid, linewidth=0.8, markersize=4, fmt='o', label='Single stars', color='royalblue')	# the error bars are too small to be visible on the plot
plt.errorbar(Prot_binaries[crit_wide_binaries], SCaHK_binaries[crit_wide_binaries], yerr=sigma_SCaHK_binaries[crit_wide_binaries], linewidth=0.8, markersize=4, fmt='o', label='Wide binaries', color='r')			# the error bars are too small to be visible on the plot
plt.errorbar(Prot_binaries[crit_close_binaries], SCaHK_binaries[crit_close_binaries], yerr=sigma_SCaHK_binaries[crit_close_binaries], linewidth=0.8, markersize=4, fmt='o', label='Close binaries', color='k')			# the error bars are too small to be visible on the plot
plt.xlim(-10, 210)
plt.xlabel(r'$P_{rot}$' + ' (days)', fontsize='x-large')
plt.ylabel(r'$S_{Ca_{II, H&K}}$', fontsize='x-large')
plt.title('Level of magnetic activity as probed by the CaII H&K\nspectral lines as a function of the rotation period of stars', fontsize='x-large')
plt.legend()
plt.savefig('./CaHK_activity_versus_rotation_single_stars_versus_binaries.pdf', format='pdf')
plt.close()


