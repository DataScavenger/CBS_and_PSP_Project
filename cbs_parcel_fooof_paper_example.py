### Create examples ->  'fixed'
#                   ->  'knee'
#                   ->  'offset'
#                   ->  'stepwise'

#_____ Import Libraries _____

#Import numpy
import numpy as np
#Import matplotlib
import matplotlib.pyplot as plt
#Import scipy
import scipy.io as sio
#Import the FOOOF object
from fooof import FOOOF
    
# Example Plot for differenct models

###############
# Fixed Model #
###############

sub = 'cbs06'
label = 'Precentral_R'

#Load the power spectrum
spectra = sio.loadmat('c:/data/parcel/power/' + sub + '/' + sub + '_parcel_pow.mat')
spectra = spectra.get('psd_for_fooof')
#for plotting we need frequencies
freqs =  np.round( spectra['freq'][0, 0][0, ] , decimals = 1 )
labels = spectra['label'][0,0] #index content with labels[0,0][0], labels[1,0][0], ...
spec = spectra['powspctrm'][0,0]
#select spectrum from 'Precentral_R'
idx = np.where(labels==label)[0][0]
spectrum = spec[idx,:]

#fit fixed model
fm_fixed = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_fixed.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[1][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_fixed.freqs,fm_fixed.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_fixed.freqs,fm_fixed.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_fixed.freqs,fm_fixed._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_fixed.gaussian_params_[0][0], fm_fixed.gaussian_params_[0][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_1], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_fixed.gaussian_params_[0][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_fixed.gaussian_params_[1][0], fm_fixed.gaussian_params_[1][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_2], fm_fixed.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_fixed.gaussian_params_[1][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
ax.set_xlim(3,48)
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/fixed_model_01.tiff')
plt.show()
plt.close()

##############
# Knee Model #
##############

sub = 'pd07'
label = 'Frontal_Inf_Tri_R'

#Load the power spectrum
spectra = sio.loadmat('c:/data/parcel/power/' + sub + '/' + sub + '_parcel_pow.mat')
spectra = spectra.get('psd_for_fooof')
#for plotting we need frequencies
freqs =  np.round( spectra['freq'][0, 0][0, ] , decimals = 1 )
labels = spectra['label'][0,0] #index content with labels[0,0][0], labels[1,0][0], ...
spec = spectra['powspctrm'][0,0]
#select spectrum from 'Frontal_Inf_Tri_R'
idx = np.where(labels==label)[0][0]
spectrum = spec[idx,:]

#1.begin with 'failed' fixed model

#fit fixed model
fm_fixed = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_fixed.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[1][0]) )
closest_freq_gauss_3 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[2][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_fixed.freqs,fm_fixed.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_fixed.freqs,fm_fixed.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_fixed.freqs,fm_fixed._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_fixed.gaussian_params_[0][0], fm_fixed.gaussian_params_[0][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_1], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_fixed.gaussian_params_[0][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_fixed.gaussian_params_[1][0], fm_fixed.gaussian_params_[1][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_2], fm_fixed.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_fixed.gaussian_params_[1][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_fixed.gaussian_params_[2][0], fm_fixed.gaussian_params_[2][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_3], fm_fixed.fooofed_spectrum_[closest_freq_gauss_3] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_fixed.gaussian_params_[2][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_3], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/knee_model_01.tiff')
plt.show()
plt.close()

#2.knee model

#fit knee model
fm_knee = FOOOF(aperiodic_mode='knee', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_knee.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[1][0]) )
closest_freq_gauss_3 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[2][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_knee.freqs,fm_knee.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_knee.freqs,fm_knee.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_knee.freqs,fm_knee._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_knee.gaussian_params_[0][0], fm_knee.gaussian_params_[0][0] ], [ fm_knee._ap_fit[closest_freq_gauss_1], fm_knee.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_knee.gaussian_params_[0][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_knee.gaussian_params_[1][0], fm_knee.gaussian_params_[1][0] ], [ fm_knee._ap_fit[closest_freq_gauss_2], fm_knee.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_knee.gaussian_params_[1][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_knee.gaussian_params_[2][0], fm_knee.gaussian_params_[2][0] ], [ fm_knee._ap_fit[closest_freq_gauss_3], fm_knee.fooofed_spectrum_[closest_freq_gauss_3] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_knee.gaussian_params_[2][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_3], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/knee_model_02.tiff')
plt.show()
plt.close()

################
# Offset Model #
################

sub = 'pd53'
label = 'Temporal_Sup_L'

#Load the power spectrum
spectra = sio.loadmat('c:/data/parcel/power/' + sub + '/' + sub + '_parcel_pow.mat')
spectra = spectra.get('psd_for_fooof')
#for plotting we need frequencies
freqs =  np.round( spectra['freq'][0, 0][0, ] , decimals = 1 )
labels = spectra['label'][0,0] #index content with labels[0,0][0], labels[1,0][0], ...
spec = spectra['powspctrm'][0,0]
#select spectrum from 'Frontal_Inf_Tri_R'
idx = np.where(labels==label)[0][0]
spectrum = spec[idx,:]

#1.begin with 'failed' fixed model

#fit fixed model
fm_fixed = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_fixed.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[1][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_fixed.freqs,fm_fixed.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_fixed.freqs,fm_fixed.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_fixed.freqs,fm_fixed._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_fixed.gaussian_params_[0][0], fm_fixed.gaussian_params_[0][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_1], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_fixed.gaussian_params_[0][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_fixed.gaussian_params_[1][0], fm_fixed.gaussian_params_[1][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_2], fm_fixed.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_fixed.gaussian_params_[1][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/offset_model_01.tiff')
plt.show()
plt.close()

#2.go on with 'failed' knee model

#fit knee model
fm_knee = FOOOF(aperiodic_mode='knee', peak_width_limits=[2,14] ,peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_knee.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[1][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_knee.freqs,fm_knee.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_knee.freqs,fm_knee.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_knee.freqs,fm_knee._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_knee.gaussian_params_[0][0], fm_knee.gaussian_params_[0][0] ], [ fm_knee._ap_fit[closest_freq_gauss_1], fm_knee.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_knee.gaussian_params_[0][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_knee.gaussian_params_[1][0], fm_knee.gaussian_params_[1][0] ], [ fm_knee._ap_fit[closest_freq_gauss_2], fm_knee.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_knee.gaussian_params_[1][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/offset_model_02.tiff')
plt.show()
plt.close()

#3.offset model

#force the offset to be on the 3Hz frequency bin
freq_range_start = 3
freq_range_end = 48
#prune freqs & spectrum to frequency range (This is necessary to force the offset in the aperiodic fit to be at 'freq_range_start' and not on the first frequency bin in freqs)
start = np.where( freqs == freq_range_start )[0][0]
stop =  np.where( freqs == freq_range_end )[0][0]
#adapt 'freqs' and 'spectrum'
freq_range = freqs[range(start,stop+1)]
offset_spectrum = spectrum[range(start,stop+1)]

#initialize the FOOOF model
fm_offset = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
fm_offset._ap_bounds=(( -np.inf,-np.inf,-np.inf ),( np.inf,np.inf,np.inf ))

#fit FOOOF model
fm_offset.fit(np.arange(1,len(freq_range)+1), offset_spectrum)            # Setting the range fom 1 to len..+1 is only to allow the fit to run without warnings / errors. Indeed, frequencies used are freq_range (in the background)
fm_offset.freqs = freq_range                                              # Set freqs
fm_offset.freq_range = [freq_range[0],freq_range[-1]]                     # Set freq_range
fm_offset.aperiodic_mode = 'offset'                                       # change to aperiodic mode "initial"

#adapt peak frequency of gaussian fits (which have been fitted on 1Hz:91Hz, instead of 3Hz:48Hz)
fm_offset.gaussian_params_[0][0] = fm_offset.gaussian_params_[0][0] * 0.5 + freq_range_start - 0.5
fm_offset.gaussian_params_[1][0] = fm_offset.gaussian_params_[1][0] * 0.5 + freq_range_start - 0.5
fm_offset.gaussian_params_[2][0] = fm_offset.gaussian_params_[2][0] * 0.5 + freq_range_start - 0.5

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_offset.freqs - fm_offset.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_offset.freqs - fm_offset.gaussian_params_[1][0]) )
closest_freq_gauss_3 = np.argmin( abs(fm_offset.freqs - fm_offset.gaussian_params_[2][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_offset.freqs,fm_offset.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_offset.freqs,fm_offset.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_offset.freqs,fm_offset._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_offset.gaussian_params_[0][0], fm_offset.gaussian_params_[0][0] ], [ fm_offset._ap_fit[closest_freq_gauss_1], fm_offset.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_offset.gaussian_params_[0][0], fm_offset.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_offset.gaussian_params_[1][0], fm_offset.gaussian_params_[1][0] ], [ fm_offset._ap_fit[closest_freq_gauss_2], fm_offset.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_offset.gaussian_params_[1][0], fm_offset.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_offset.gaussian_params_[2][0], fm_offset.gaussian_params_[2][0] ], [ fm_offset._ap_fit[closest_freq_gauss_3], fm_offset.fooofed_spectrum_[closest_freq_gauss_3] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_offset.gaussian_params_[2][0], fm_offset.fooofed_spectrum_[closest_freq_gauss_3], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/offset_model_03.tiff')
plt.show()
plt.close()

##################
# Stepwise Model #
##################

sub = 'pd26'
label = 'Frontal_Inf_Orb_L'

#Load the power spectrum
spectra = sio.loadmat('c:/data/parcel/power/' + sub + '/' + sub + '_parcel_pow.mat')
spectra = spectra.get('psd_for_fooof')
#for plotting we need frequencies
freqs =  np.round( spectra['freq'][0, 0][0, ] , decimals = 1 )
labels = spectra['label'][0,0] #index content with labels[0,0][0], labels[1,0][0], ...
spec = spectra['powspctrm'][0,0]
#select spectrum from 'Frontal_Inf_Tri_R'
idx = np.where(labels==label)[0][0]
spectrum = spec[idx,:]

#1.begin with 'failed' fixed model

#fit fixed model
fm_fixed = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_fixed.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_fixed.freqs - fm_fixed.gaussian_params_[0][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_fixed.freqs,fm_fixed.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_fixed.freqs,fm_fixed.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_fixed.freqs,fm_fixed._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_fixed.gaussian_params_[0][0], fm_fixed.gaussian_params_[0][0] ], [ fm_fixed._ap_fit[closest_freq_gauss_1], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_fixed.gaussian_params_[0][0], fm_fixed.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/stepwise_model_01.tiff')
plt.show()
plt.close()

#2.knee model

#fit knee model
fm_knee = FOOOF(aperiodic_mode='knee', peak_width_limits=[2,14] ,peak_threshold=2, max_n_peaks=4)
#fit FOOOF model
fm_knee.fit(freqs, spectrum,[3,48])

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_knee.freqs - fm_knee.gaussian_params_[1][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_knee.freqs,fm_knee.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_knee.freqs,fm_knee.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_knee.freqs,fm_knee._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_knee.gaussian_params_[0][0], fm_knee.gaussian_params_[0][0] ], [ fm_knee._ap_fit[closest_freq_gauss_1], fm_knee.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_knee.gaussian_params_[0][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_knee.gaussian_params_[1][0], fm_knee.gaussian_params_[1][0] ], [ fm_knee._ap_fit[closest_freq_gauss_2], fm_knee.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_knee.gaussian_params_[1][0], fm_knee.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/stepwise_model_02.tiff')
plt.show()
plt.close()

#3.offset model

#force the offset to be on the 3Hz frequency bin
freq_range_start = 3
freq_range_end = 48
#prune freqs & spectrum to frequency range (This is necessary to force the offset in the aperiodic fit to be at 'freq_range_start' and not on the first frequency bin in freqs)
start = np.where( freqs == freq_range_start )[0][0]
stop =  np.where( freqs == freq_range_end )[0][0]
#adapt 'freqs' and 'spectrum'
freq_range = freqs[range(start,stop+1)]
offset_spectrum = spectrum[range(start,stop+1)]

#initialize the FOOOF model
fm_offset = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14], peak_threshold=2, max_n_peaks=4)
fm_offset._ap_bounds=(( -np.inf,-np.inf,-np.inf ),( np.inf,np.inf,np.inf ))

#fit FOOOF model
fm_offset.fit(np.arange(1,len(freq_range)+1), offset_spectrum)            # Setting the range fom 1 to len..+1 is only to allow the fit to run without warnings / errors. Indeed, frequencies used are freq_range (in the background)
fm_offset.freqs = freq_range                                               # Set freqs
fm_offset.freq_range = [freq_range[0],freq_range[-1]]                      # Set freq_range
fm_offset.aperiodic_mode = 'offset'                                        # change to aperiodic mode "initial"

#adapt peak frequency of gaussian fits
fm_offset.gaussian_params_[0][0] = fm_offset.gaussian_params_[0][0] * 0.5 + freq_range_start - 0.5
fm_offset.gaussian_params_[1][0] = fm_offset.gaussian_params_[1][0] * 0.5 + freq_range_start - 0.5

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
closest_freq_gauss_1 = np.argmin( abs(fm_offset.freqs - fm_offset.gaussian_params_[0][0]) )
closest_freq_gauss_2 = np.argmin( abs(fm_offset.freqs - fm_offset.gaussian_params_[1][0]) )
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(fm_offset.freqs,fm_offset.power_spectrum,color='black',linewidth=4,label='Original Spectrum')
ax.plot(fm_offset.freqs,fm_offset.fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
ax.plot(fm_offset.freqs,fm_offset._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
ax.plot( [ fm_offset.gaussian_params_[0][0], fm_offset.gaussian_params_[0][0] ], [ fm_offset._ap_fit[closest_freq_gauss_1], fm_offset.fooofed_spectrum_[closest_freq_gauss_1] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
ax.scatter( fm_offset.gaussian_params_[0][0], fm_offset.fooofed_spectrum_[closest_freq_gauss_1], s=350, marker='o', color='springgreen' )
ax.plot( [ fm_offset.gaussian_params_[1][0], fm_offset.gaussian_params_[1][0] ], [ fm_offset._ap_fit[closest_freq_gauss_2], fm_offset.fooofed_spectrum_[closest_freq_gauss_2] - 0.015], linewidth=4, color='springgreen')
ax.scatter( fm_offset.gaussian_params_[1][0], fm_offset.fooofed_spectrum_[closest_freq_gauss_2], s=350, marker='o', color='springgreen' )
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/stepwise_model_03.tiff')
plt.show()
plt.close()

#4.stepwise model

FooofModels = []
steps = [3.,15,48.]

# Set initial offset for model fit [as last power of the last freqbin of the previous fit, if there has been one]
offset = np.log10(spectrum[np.where(freqs == steps[0])[0][0]])

#fit (#steps - 1) fooof models
for k in range(len(steps)-1):
    #get borders for kth interval
    von = steps[k]
    bis = steps[k+1]
    #identify borders in 'freqs' and set frequency range upon which to fit the model
    start = np.where( freqs == von )[0][0]
    stop = np.where( freqs == bis )[0][0]
    freq_range = freqs[range(start,stop + 1)]
    #get kth intervall-spectrum to fit
    tmp_spectrum = spectrum[range(start,stop + 1)]
    
    #fit Model
    tmp = FOOOF(aperiodic_mode='fixed', peak_width_limits=[2,14],peak_threshold=2, max_n_peaks=4)
    tmp._ap_bounds= ((-np.inf,-np.inf,-np.inf),(np.inf,np.inf,np.inf))
    
    #store single Fooof-Models in list
    FooofModels.append(tmp)
    #fit fixed model
    tmp.fit(np.arange(1,len(freq_range)+1), tmp_spectrum)   # Setting the range fom 1 to len..+1 is only to allow the fit to run without warnings / errors. Indeed, frequencies used are freq_range (in the background)
    tmp.freqs = freq_range                                  # set freqs 
    tmp.freq_res = freq_range[1] - freq_range[0]            # freq resolution
    tmp.freq_range = [start,stop]                           # set freq_range
    
    #set new offset to be the last
    offset = tmp.power_spectrum[-1]
    
    #clear tmp
    del(tmp)

#create figure
fig, ax = plt.subplots(figsize=[22.5, 12])
plt.rcParams['font.size'] = '45'
plt.rcParams['font.weight'] = 'bold'
ax.plot(freqs,np.log10(spectrum), color='black',linewidth=4,label='Original Spectrum')
for k in range(len(FooofModels)):
    for g in range(len(FooofModels[k].gaussian_params_)):
        #adapt gaussian
        FooofModels[k].gaussian_params_[g][0] = FooofModels[k].gaussian_params_[g][0] * 0.5 + steps[k] - 0.5
        closest_freq_gauss_g = np.argmin( abs( FooofModels[k].freqs - FooofModels[k].gaussian_params_[g][0]) )
        #plot gth gaussian
        ax.plot( [ FooofModels[k].gaussian_params_[g][0], FooofModels[k].gaussian_params_[g][0] ], [ FooofModels[k]._ap_fit[closest_freq_gauss_g], FooofModels[k].fooofed_spectrum_[closest_freq_gauss_g] - 0.015], linewidth=4, color='springgreen')
        ax.scatter( FooofModels[k].gaussian_params_[g][0], FooofModels[k].fooofed_spectrum_[closest_freq_gauss_g], s=350, marker='o', color='springgreen' )
    #go on with remaining plot
    if k == len(FooofModels)-1:
        ax.plot(FooofModels[k].freqs,  FooofModels[k].fooofed_spectrum_,color='indianred',linewidth=4,label='Fooof Model Fit')
        ax.plot(FooofModels[k].freqs,  FooofModels[k]._ap_fit,'--',color='royalblue',linewidth=4,label='Aperiodic Fit')
        #plot gth gaussian
        ax.plot( [ FooofModels[k].gaussian_params_[g][0], FooofModels[k].gaussian_params_[g][0] ], [ FooofModels[k]._ap_fit[closest_freq_gauss_g], FooofModels[k].fooofed_spectrum_[closest_freq_gauss_g] - 0.015], linewidth=4, color='springgreen', label='Peak Fit')
        ax.scatter( FooofModels[k].gaussian_params_[g][0], FooofModels[k].fooofed_spectrum_[closest_freq_gauss_g], s=350, marker='o', color='springgreen' )
    else:
        ax.plot(FooofModels[k].freqs,  FooofModels[k].fooofed_spectrum_,color='indianred',linewidth=4)
        ax.plot(FooofModels[k].freqs,  FooofModels[k]._ap_fit,'--',color='royalblue',linewidth=4)
ax.set_xlim(3,48)
ax.set_ylabel('log10(Power)',fontweight='bold')
ax.set_xlabel('Frequency',fontweight='bold')
leg = ax.legend()
fig.savefig('c:/data/fooof_flexible/stepwise_model_04.tiff')
plt.show()
plt.close()



