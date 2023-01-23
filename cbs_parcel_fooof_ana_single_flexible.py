"""
Script for splitting PSD into periodic & aperiodic components & peak estimates

"""

#_____ Import Libraries _____

#Import numpy
import numpy as np
#Import matplotlib
import matplotlib.pyplot as plt
#Import scipy
import scipy.io as sio
#Import os
import os
#Import the FOOOF object
from fooof import FOOOF

#################################################################
############# __________F_U_N_C_T_I_O_N_S__________ ############# -> Needed for the Script. Function to identify Troughs
#################################################################

def my_ismember(list1,list2):
    
    """ This function works like the matlab function 'ismember()'. This said,
        it checks which elements in list1 are also present in list2. If an element is
        present it is set to 'True' at its position/index in the array list1, elsewise
        it is set to 'False'. """
        
    list1= list(list1)
    list2 = list(list2)
    boolean = []
    for i in range(len(list1)):
        if list1[i] in list2:
            boolean.extend([True])
        else:
            boolean.extend([False])
    return boolean
    
    
def my_r_squared(FooofModel):
    
    """ This function takes a single Fooof-Model and (eventually) multiple Fooof-Models forwarded in a list
        and calculates the R² as it is described in 'Introduction into Statistical Learnung' p.70 with the use
        of Residual Sum of Squares (RSS) and Total Sum of Squares (TSS). This fits to the values which are given
        by a single Fooof-Model with e.g. fm.r_squared_ """
    #check if it is already a list as the input
    if isinstance(FooofModel,list) == False:
        L = [FooofModel]
    else:
        L = FooofModel
    #get Fooof-Spectrum and Original Power-Spectrum
    FS = []
    PS = []
    for k in range(len(L)):
        if k == len(L)-1:
            FS.extend( L[k].fooofed_spectrum_[ range(len(L[k].fooofed_spectrum_)) ] )
            PS.extend( L[k].power_spectrum[ range(len(L[k].power_spectrum)) ] )
        else:
            FS.extend( L[k].fooofed_spectrum_[ range(len(L[k].fooofed_spectrum_)-1) ] )
            PS.extend( L[k].power_spectrum[ range(len(L[k].power_spectrum)-1) ] )
    #calculate RSS and TSS
    RSS = sum( np.square( np.array(PS) - np.array(FS) ) )
    TSS = sum( np.square( np.array(PS) - (np.sum(PS)/(len(PS))) ) )
    #get R²
    R_square = 1 - (RSS/TSS)
    
    return R_square
    
    
def my_troughs(spectrum_flat, spectrum_freqs, threshold, trough_intervall, show):
    
    """ This function returns frequency bins of the input data that subseed the given threshold (troughs within intervall, cliff).
        The input needed is the 'flat spectrum' of the output of a fooof-model fitted to a single spectrum
        and the 'frequencies' used as a range to fit this model.
        The function outputs a tupel:
        The first output array constitutes frequency bins within 'frequencies', where the first element is the
        initial frequency bin used for fitting the model, followed by the frequency bins of eventual troughs
        and finally the last frequency bin used for fitting the model. E.g.
        array = [ freq_range[0], troughs ,freq_range[last] ].
        The second output variable in the tupel is a logical indicating if a 'Cliff' e.g. the initial power value
        at the first frequency bin is further than the threshold away from the aperiodic model fit."""
    
    # Get Zero-Crossings Indices & corresponsing frequency bins
    zero_cross = np.where( spectrum_flat < 0 )[0]
    zero_cross_freq = spectrum_freqs[ spectrum_flat < 0 ]
    # Visualize potential 'Issues'
    if show == True:
        plt.figure()
        plt.plot( spectrum_freqs, spectrum_flat )
        plt.scatter( spectrum_freqs, spectrum_flat )
        plt.scatter( zero_cross_freq, np.zeros(np.shape(zero_cross)) )
        plt.plot( spectrum_freqs, np.zeros(np.shape(spectrum_flat)) )
        plt.scatter( zero_cross_freq, spectrum_flat[zero_cross] )
        plt.show()
    
    # Find number of Zero Crossing Intervalls [new negative range begins with a non-1 number]
    num_intervalls = np.concatenate( (np.array([1]),np.diff(zero_cross)), axis = 0) != 1
    if np.size(zero_cross) != 0:
        # Beginning of Intervall [Start Index]
        StartIndex = np.concatenate( ([zero_cross[0]], zero_cross[num_intervalls]) )
        # End of Intervall [End Index]
        EndIndex =  np.concatenate( (zero_cross[0:-1][ num_intervalls[1:len(num_intervalls)] ], np.array([zero_cross[-1]]) ) )
    else:
        StartIndex = np.array([])
        EndIndex = np.array([])
    
    # Troughs of each Intervall
    Troughs = []
    for i in range(len(StartIndex)):
        # Indices Intervall
        intervall = np.array(range(StartIndex[i],EndIndex[i]+1))
        # Get Troughs (:= minima within the intervalls & subseeding threshold)
        if abs(min(spectrum_flat[intervall])) >= threshold[0]:
            Troughs.extend( [intervall[ np.argmin( spectrum_flat[intervall] ) ]] )
    
    # Add initial freqbin and final freqbin
    Troughs.extend([0,-1])
    Troughs = list(set(Troughs))
    Troughs = spectrum_freqs[Troughs]
    Troughs.sort()                                                      # => Trough now contains offset at [0] and last frequency bin to fit ['end']
    # apply subseed_range
    delete_Trough = ( Troughs < trough_intervall[0] ) |	 ( Troughs > trough_intervall[-1] )
    delete_Trough[ [0,-1], ] = False                                    # keep initial freqbin and last freqbin
    # clear Troughs
    Troughs = np.delete(Troughs,delete_Trough)
    
    # Check for Cliff
    if abs( spectrum_flat[0] ) > threshold[1]:
        Cliff = True
    else:
        Cliff = False
    
    return Troughs, Cliff
    
    
def my_thresh(spectrum_flat, spectrum_freqs, threshold, thresh_intervall, show):
    
    """ This function returns frequency bins of the input data that sub/exceed the given threshold (troughs within intervall, cliff).
        The input needed is the 'flat spectrum' of the output of a fooof-model fitted to a single spectrum
        and the 'frequencies' used as a range to fit this model.
        The function outputs a tupel:
        The first output array constitutes frequency bins within 'frequencies', where the first element is the
        initial frequency bin used for fitting the model, followed by the frequency bins of eventual troughs
        and finally the last frequency bin used for fitting the model. E.g.
        array = [ freq_range[0], troughs ,freq_range[last] ].
        The second output variable in the tupel is a logical indicating if a 'Cliff' e.g. the initial power value
        at the first frequency bin is further than the threshold away from the aperiodic model fit."""
    
    # Get sub/exceeding Threshold Indices & corresponding frequency bins
    thresh_cross = np.where( spectrum_flat < -threshold[0] )[0]
    thresh_cross_freq = spectrum_freqs[ spectrum_flat < -threshold[0] ]
    # Restrict Intervall to thresh_intervall
    step = spectrum_freqs[1] - spectrum_freqs[0]
    keepthese = my_ismember(thresh_cross_freq, np.arange(thresh_intervall[0], thresh_intervall[-1]+step, step))#range(thresh_intervall[0],thresh_intervall[-1]+1))
    thresh_cross = thresh_cross[keepthese]
    thresh_cross_freq = thresh_cross_freq[keepthese]
    
    # Visualize potential 'Issues'
    if show == True:
        plt.figure()
        #flat spectrum
        plt.plot( spectrum_freqs, spectrum_flat )
        plt.scatter( spectrum_freqs, spectrum_flat )
        #freqbins (sub/exceeding)
        plt.plot( spectrum_freqs, np.zeros(np.shape(spectrum_flat)) )
        plt.scatter( thresh_cross_freq, np.zeros(np.shape(thresh_cross)) )
        #threshold lines
        plt.plot( spectrum_freqs, np.repeat(-threshold[0],len(spectrum_freqs)),linestyle='dashed',color='grey')
        #sub/exceed spectrum values
        plt.scatter( thresh_cross_freq, spectrum_flat[thresh_cross] )
        plt.show()
    
    # Find number of Zero Crossing Intervalls [new negative range begins with a non-1 number]
    num_intervalls = np.concatenate( (np.array([1]),np.diff(thresh_cross)), axis = 0) != 1
    if np.size(thresh_cross) != 0:
        # Beginning of Intervall [Start Index]
        StartIndex = np.concatenate( ([thresh_cross[0]], thresh_cross[num_intervalls]) )
        # End of Intervall [End Index]
        EndIndex =  np.concatenate( (thresh_cross[0:-1][ num_intervalls[1:len(num_intervalls)] ], np.array([thresh_cross[-1]]) ) )
    else:
        StartIndex = np.array([])
        EndIndex = np.array([])
    
    # Troughs of each Intervall
    Troughs = []
    for i in range(len(StartIndex)):
        # Indices Intervall
        intervall = np.array(range(StartIndex[i],EndIndex[i]+1))
        # Get Troughs (:= minima within the intervalls & subseeding threshold)
        if abs(min(spectrum_flat[intervall])) >= threshold[0]:
            Troughs.extend( [intervall[ np.argmin( spectrum_flat[intervall] ) ]] )
    
    # Add initial freqbin and final freqbin
    Troughs.extend([0,-1])
    Troughs = list(set(Troughs))
    Troughs = spectrum_freqs[Troughs]
    Troughs.sort()                                                      # => Trough now contains offset at [0] and last frequency bin to fit ['end']
    
    # Check for Cliff
    if abs( spectrum_flat[0] ) > threshold[1]:
        Cliff = True
    else:
        Cliff = False
    
    return Troughs, Cliff
    
    
def my_read_txt(path):
    
    text_data = list()
    current_file = os.path.abspath(os.path.join(path))
    
    if os.path.exists(current_file):
        open_file = open(current_file, 'r', encoding="latin-1")
        text_data = open_file.read().split('\n')
        text_data = list(filter(None, text_data))
    return text_data
    
    
def my_search(response, subject, area):
    
    #store parameters
    SubList = np.empty( ( len(response),1 ) )
    ParList = np.empty( ( len(response),1 ) )
    SubList[:] = np.NaN
    ParList[:] = np.NaN
    
    for k in range(len(response)):
        
        r = response[k].split('\t')
        
        if subject in r:
            SubList[k] = r.index(subject)
        if area in r:
            ParList[k] = r.index(area)
        
    return np.where( np.isnan( SubList + ParList ) == False )[0].astype(int)
    
    
############################################################
############# __________S_C_R_I_P_T___________ #############
############################################################

### Fitting FOOOF Models

#Mainpath
mpath = 'c:/data'

#Get Txt-File Infos [fooof_flexible_settings]
txt = my_read_txt( mpath + '/fooof_flexible_settings.txt' )

#Loop variables
subjects = os.listdir(mpath + '/parcel/power')

UsedModels = []

#subject loop
for s in np.arange(len(subjects)):
    
    sub = subjects[s]   #Load subject
    
    if not os.path.exists(mpath + '/fooof_flexible/' + sub + '/images/parcel'):
        os.makedirs(mpath + '/fooof_flexible/' + sub + '/images/parcel')          #Put images
    if not os.path.exists(mpath + '/fooof_flexible/' + sub + '/results_fooof/parcel'):
        os.makedirs(mpath + '/fooof_flexible/' + sub + '/results_fooof/parcel')   #Put results
    
    #Load the power spectrum
    spectra = sio.loadmat(mpath + '/parcel/power/' + sub + '/' + sub + '_parcel_pow.mat')
    spectra = spectra.get('psd_for_fooof')
    
    # for plotting we need frequencies
    freqs =  np.round( spectra['freq'][0, 0][0, ] , decimals = 1 )
    labels = spectra['label'][0,0]
    spec = spectra['powspctrm'][0,0]
    
    #Let's see how the Power Spectra look BEFORE FOOOF modeling
    plt.figure()                                          #new figure
    plt.plot(freqs, np.transpose(np.log10(spec)))         #plot original spectra (log10)
    plt.ylabel('Power (log10)')
    plt.xlabel('Frequencies')
    plt.xlim(freqs[0],freqs[-1])
    plt.title(sub + 'Parcel Power Spectra')
    #Save image
    plt.savefig(mpath + '/fooof_flexible/' + sub + '/images/parcel/' + sub + '_parcel_spectra.png')
    #Close figure
    plt.close()
    
    #Initiaite dict to store information
    all_dict = {}
    
    #Get R² of all possible model fits
    Rsq = []
    
    #spectrum loop
    for i in range(len(spec)):
        
        spectrum = spec[i, :]    #Get Spectrum
        
        ####################################################################################################
        #_______________________________Model Specification [ either I. or II.]____________________________#
        #                                                                                                  #
        # I. Spectrum is recognized in 'txt' -> Use noted information for proper analysis                  #
        # II.Spectrum is not recognized in 'txt' -> Explore the spectrum with an algorithmic way [troughs] #
        #                                                                                                  #
        ####################################################################################################
        
        # either go 'fixed' | either go 'knee' | either go 'initial' | either go 'stepwise'
        found = my_search(txt, sub, labels[i][0][0])
        
        #read content
        if len(found) == 1:
            r = txt[found[0]]
        elif len(found) == 0:
            r = 'algo'
            print('no such result found in .txt. Use "algorithm" instead')
        else:
            print('multiple results found')
        
        #get content of 'r'
        r = r.split('\t')
        
        #read out parameter information [important: third element must be odered as freq range ([]) -> ap_bounds_ ([],[],[]) -> (eventually) peak_width_params ([])]
        if len(r) == 3:
            p = r[2].split(' ')
            if p[0] == 'fixed' or p[0] == 'knee' or p[0] == 'initial' or p[0] == 'stepwise':
                model = p[0]
                params = []
                for n in range(1,len(p)):
                    params.append( eval(p[n]) )
            else:
                print('no model was specified')
        elif len(r) == 1 and r[0] == 'algo':
            model = 'algo'
            params = [[3, 48], [-np.inf, np.inf], [-np.inf, 20], [-np.inf, np.inf], [2,14]] #default parameter settings
        else:
            print('something is wrong in r: ' + labels[i][0][0] + ' with the .txt')
        
        ####################################################################################################
        #_______________________________________Model specified____________________________________________#
        ####################################################################################################
        
        
        ####################################################################################################
        #____________________________________________Model Fit_____________________________________________#
        #                                                                                                  #
        # The following code uses prespecified model fitting as given by the .txt                          #
        # with keywords: ['fixed', 'knee', 'initial' or 'stepwise']                                        #
        # OR it uses 'algo', an algorithmic procedure where it begins with 'fixed' -> 'knee' -> 'initial'  #
        # -> 'stepwise' and succesive evaluation of each step via the sub-Zero-Line Troughs                #
        #                                                                                                  #
        ####################################################################################################
        
        if model == 'algo' or model == 'fixed':
                       
            print('algo or fixed')
                        
            #########################################
            # 1.1 'fixed' Model fit [algo or fixed] #
            #########################################
            
            # Set the frequency range upon which to fit FOOOF [also used for possible later fits]
            freq_range_start = params[0][0]
            freq_range_end = params[0][1]
            # _ap_bounds
            up = []
            low = []
            for k in range(1,len(params)-1): low.append( params[k][0] ), up.append( params[k][1] )
                        
            # Initialize the FOOOF model
            fm_fixed = FOOOF(aperiodic_mode='fixed', peak_width_limits=params[-1],
                             peak_threshold=2, max_n_peaks=4)
            
            # Fit FOOOF model
            fm_fixed.fit(freqs, spectrum,[freq_range_start,freq_range_end])
            
            fm_fixed.plot()
            
            plt.scatter( [fm_fixed.freqs[np.where(fm_fixed.freqs == 4)[0][0]], fm_fixed.freqs[np.where(fm_fixed.freqs == 30)[0][0]]],
                         [fm_fixed.power_spectrum[np.where(fm_fixed.freqs == 4)[0][0]], fm_fixed.power_spectrum[np.where(fm_fixed.freqs == 30)[0][0]]], s = 500 ,c = 'black', marker = '*')
            plt.title('1. Fixed Model Fit')
            
            plt.savefig(mpath + '/fooof_flexible/' + sub + '/images/parcel/' + sub + labels[i][0][0] + '_fixed_model2.png')
            plt.close()
            
            if model == 'fixed':
                # ... save it
                print('save fixed model')
                # ... overwrite 'model', that the other 'if-else' statements are not executed
                model = fm_fixed
                #protocoll in UsedModels
                UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'fixed')
            else:
                
                ########################################
                # 1.2 get troughs & evaluate model fit #
                ########################################
                
                spectrum_flat = fm_fixed._spectrum_flat
                spectrum_freqs = fm_fixed.freqs
                # settings
                threshold = [0.05, 0.1] #threshold for troughs, threshold for cliff
                trough_intervall = np.array([7,25])
                show = False
                
                # get troughs (array of frequency bins) and cliff (True/False)
                Troughs_fixed = my_thresh(spectrum_flat, spectrum_freqs, threshold, trough_intervall, show)
                
                # get R² for fixed model
                Rsq.extend( [my_r_squared(fm_fixed)] )
                
                if ( len( Troughs_fixed[0] ) <= 2 ) & ( Troughs_fixed[1] == False ):
                    #... save it
                    print('save fixed model')
                    # ... overwrite 'model', that the other 'if-else' statements are not executed
                    model = fm_fixed
                    #protocoll in UsedModels
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'fixed')
                else:
                    print('go on with knee')
                    
        if model == 'algo' or model == 'knee':
            
            print('algo or knee')
            
            #######################################
            # 2.1 'knee' Model fit [algo or knee] #
            #######################################
            
            # Set the frequency range upon which to fit FOOOF [also used for possible later fits]
            freq_range_start = params[0][0]
            freq_range_end = params[0][1]
            # _ap_bounds
            up = []
            low = []
            for k in range(1,len(params)-1): low.append( params[k][0] ), up.append( params[k][1] )
            
            # Initialize the FOOOF model
            fm_knee = FOOOF(aperiodic_mode='knee', peak_width_limits=params[-1],
                       peak_threshold=2, max_n_peaks=4)
            fm_knee._ap_bounds = ( tuple(low),tuple(up) )
            
            # Fit FOOOF model
            fm_knee.fit(freqs,spectrum,[freq_range_start,freq_range_end])
            fm_knee.plot()
            plt.scatter( [fm_knee.freqs[np.where(fm_knee.freqs == 4)[0][0]], fm_knee.freqs[np.where(fm_knee.freqs == 30)[0][0]]],
                         [fm_knee.power_spectrum[np.where(fm_knee.freqs == 4)[0][0]], fm_knee.power_spectrum[np.where(fm_knee.freqs == 30)[0][0]]], s = 500 ,c = 'black', marker = '*')
            plt.title('2. Knee Model Fit')
            plt.savefig(mpath + '/fooof_flexible/' + sub + '/images/parcel/' + sub + labels[i][0][0] + '_knee_model2.png')
            plt.close()
            
            if model == 'knee':
                # ...  save it
                print('save knee model')
                # ... overwrite 'model', that the other 'if-else' statements are not executed
                model = fm_knee
                #protocoll in UsedModels
                UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'knee')
            else:
                
                ########################################
                # 2.2 get troughs & evaluate model fit #
                ########################################
                
                spectrum_flat = fm_knee._spectrum_flat
                spectrum_freqs = fm_knee.freqs
                # settings
                threshold = [0.05 ,0.1] #threshold for troughs, threshold for cliff
                trough_intervall = np.array([7,25])
                show = False
                
                # get troughs (array of frequency bins) and cliff (True/False)
                Troughs_knee = my_thresh(spectrum_flat, spectrum_freqs, threshold, trough_intervall, show)
                
                # get R² for knee model
                Rsq.extend( [my_r_squared(fm_knee)] )
                
                # check 'my_trough' results
                if ( len( Troughs_knee[0] ) <= 2 ) & ( Troughs_knee[1] == False ):
                    # save the stuff as before ...
                    print('save knee model')
                    model = fm_knee
                    #protocoll in UsedModels
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'knee')
                else:
                    print('go on with initial')
                
        if model == 'algo' or model == 'initial':
            
            print('algo or initial')
            
            #############################################
            # 3.1 'initial' Model fit [algo or initial] #
            #############################################
            
            # settings
            threshold = [0.05, 0.1] #threshold for troughs, threshold for cliff
            trough_intervall = np.array([7,25])
            # Set the frequency range upon which to fit FOOOF [also used for possible later fits]
            freq_range_start = params[0][0]
            freq_range_end = params[0][1]
            # Prune freqs & spectrum to frequency range (This is necessary to force the offset in the aperiodic fit to be at the 'freq_range_start' and not at '1' on the x-axis, which would be default)
            start = np.where( freqs == freq_range_start )[0][0]
            stop =  np.where( freqs == freq_range_end )[0][0]
            # get freq range and spectrum
            freq_range = freqs[range(start,stop+1)]
            initial_spectrum = spectrum[range(start,stop+1)]
            # offset
            if model == 'algo':
                lower_bound = np.log10(initial_spectrum[0]) - threshold[0]
                upper_bound = np.log10(initial_spectrum[0]) + threshold[0]
            else:
                lower_bound = params[1][0]
                upper_bound = params[1][1]
            # Initialize the FOOOF model [with bound offset]
            fm_initial = FOOOF(aperiodic_mode='fixed', peak_width_limits=params[-1],
                               peak_threshold=2, max_n_peaks=4)
            fm_initial._ap_bounds=(( lower_bound,-np.inf,-np.inf ),( upper_bound, np.inf, np.inf ))
            
            # Some plotting of the offset intervall
            plt.figure(0)
            plt.plot( freq_range, np.log10(initial_spectrum) )
            plt.scatter( freq_range[0],np.log10(initial_spectrum[0]) )
            plt.scatter( [ freq_range[0] , freq_range[0] ], [lower_bound,upper_bound] )
            plt.scatter( freq_range[0], fm_initial.aperiodic_params_[0] )
            plt.close()
            
            try:
                # Fit FOOOF model
                fm_initial.fit(np.arange(1,len(freq_range)+1), initial_spectrum)            # Setting the range fom 1 to len..+1 is only to allow the fit to run without warnings / errors. Indeed, frequencies used are freq_range (in the background)
                fm_initial.freqs = freq_range                                               # Set freqs
                fm_initial.freq_range = [freq_range[0],freq_range[-1]]                      # Set freq_range
                fm_initial.aperiodic_mode = 'initial'                                       # change to aperiodic mode "initial"
                # Some Plotting
                fm_initial.plot()                                                               # Fooof Model Fit
                plt.plot(freq_range,np.repeat(fm_initial.aperiodic_params_[0],len(freq_range))) # Offset
                plt.plot(freqs,np.log10(spectrum))                                              # Original Spectrum over all frequencies
                plt.title('3. Initial StepWise Fixed Model Fit')
                plt.savefig(mpath + '/fooof_flexible/' + sub + '/images/parcel/' + sub + labels[i][0][0] + '_StepWise_initial_fixed_model2.png')
                plt.close()
                
                if model == 'initial':
                    # ... save it
                    print('save initial model')
                    # ... overwrite 'model', that the other 'if-else' statements are not executed
                    model = fm_initial
                    #protocoll in UsedModels
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'initial')
                else:
                    
                    ########################################
                    # 3.2 get troughs & evaluate model fit #
                    ########################################
                    
                    # spectral properties
                    spectrum_flat = fm_initial._spectrum_flat
                    spectrum_freqs = fm_initial.freqs
                    # settings
                    threshold = [0.05 ,0.1] #threshold for troughs, threshold for cliff
                    trough_intervall = np.array([7,25])
                    show = False
                    
                    # ... or compute troughs ... and evalute model
                    Troughs_initial = my_troughs(spectrum_flat, spectrum_freqs, threshold, trough_intervall, show)
                    
                    # get R² for knee model
                    Rsq.extend( [my_r_squared(fm_initial)] )
                    
                    if ( len( Troughs_initial[0] ) <= 2 ) & ( Troughs_initial[1] == False ):
                        # ... save the initial model
                        print('save initial model')
                        model = fm_initial
                        #protocoll in UsedModels
                        UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'initial')
                    else:
                        print('go on with stepwise')
                    
            except:
                #initial fit failed -> no troughs for stepwise fit
                print('Error in Initial-Fooof',str(i),'-> Keep Fixed / Knee (R² for decision) & pass to next iteration')
                # either or it must be ... if in doubt choose the one with higher R².
                if Rsq[0] > Rsq[1]:
                    model = fm_fixed
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'fixed')
                else:
                    model = fm_knee
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'knee')
                    
        if model == 'algo' or model == 'stepwise':
            
            print('algo or stepwise')
            
            ###############################################
            # 4.1 'Stepwise' Model fit [algo or StepWise] #
            ###############################################
            
            try:
                
                # Store single Fooof-Models in list
                FooofModels = []
                # Extract Troughs only
                if model == 'stepwise':
                    Troughs = params[1]
                elif model == 'algo':
                    Troughs = Troughs_initial[0]
                else:
                    UsedModels.append('Failed StepWise-Fit')
                
                # Set initial offset for model fit [as last power of the last freqbin of the previous fit, if there has been one]
                offset = np.log10(spectrum[np.where(freqs == Troughs[0])[0][0]])
                
                threshold = [0.05, 0.1] #threshold for for the stepwise fitting procedure
                
                # Iteratively fit From Offset -> Trough -> End => Check if needs to be repeated
                for k in range(len(Troughs)-1):
                    # Set the frequency range upon which to fit FOOOF
                    von = Troughs[k]
                    bis = Troughs[k+1]
                    # start & stop indices & freq_range
                    start = np.where( freqs == von )[0][0]
                    stop = np.where( freqs == bis )[0][0]
                    freq_range = freqs[range(start,stop + 1)]
                    # Get 'Intervall'-Spectrum to fit
                    tmp_spectrum = spectrum[range(start,stop + 1)]
                    
                    # Fit Model
                    tmp = FOOOF(aperiodic_mode='fixed', peak_width_limits=params[-1],
                                peak_threshold=2, max_n_peaks=4)
                    tmp._ap_bounds= (( offset - threshold[0]/2,-np.inf,-np.inf),
                                     ( offset + threshold[0]/2, np.inf, np.inf))
                    
                    # Store single Fooof-Models in list
                    FooofModels.append(tmp)
                    # Fit FOOOF model
                    tmp.fit(np.arange(1,len(freq_range)+1), tmp_spectrum)   # Setting the range fom 1 to len..+1 is only to allow the fit to run without warnings / errors. Indeed, frequencies used are freq_range (in the background)
                    tmp.freqs = freq_range                                  # set freqs 
                    tmp.freq_res = freq_range[1] - freq_range[0]            # freq resolution
                    tmp.freq_range = [start,stop]                           # set freq_range
                    tmp.aperiodic_mode = 'stepwise'
                    tmp.plot()
                    plt.scatter( tmp.freqs[-1],tmp._ap_fit[-1] )
                    plt.show()
                    
                    # set new offset
                    offset = tmp.power_spectrum[-1]
                    
                    # clear tmp
                    del(tmp)
                    
                # investiage result of Step Wise fitting
                plt.plot(freqs,np.log10(spectrum), color='black')
                for k in range(len(FooofModels)):
                    #get model
                    m = FooofModels[k]
                    #plot properties
                    plt.plot(m.freqs, m.fooofed_spectrum_, linestyle='dashed', color='tomato')
                    plt.plot(m.freqs, m._ap_fit, linestyle='dashed', color='cornflowerblue')
                    plt.title('3. Complete StepWise Fixed Model Fit')
                plt.savefig(mpath + '/fooof_flexible/' + sub + '/images/parcel/' + sub + labels[i][0][0] + '_StepWise_complete_fixed_model2.png')
                plt.close()
                print('save the StepWise model and stuff')
                
                # save the StepWise Model
                model = FooofModels
                
                #protocoll in UsedModels
                UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'StepWise')
                
                # get R² for StepWise Model
                Rsq.extend( [my_r_squared(FooofModels)] )                                       
                
            except:
                #StepWise fit failed
                print('Error in StepWise-Fooof',str(i),'-> Keep Fixed / Knee (R² for decision) & pass to next iteration')
                # either or it must be ... if in doubt choose the one with higher R².
                if Rsq[0] > Rsq[1]:
                    model = fm_fixed
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'fixed')
                else:
                    model = fm_knee
                    UsedModels.append(sub + ' ' + labels[i][0][0] + ' ' + 'knee')      
                
        plt.close('all')
        
        ####################################################################################################
        #_________________________________________END Model Fit____________________________________________#
        #                                                                                                  #
        # The variable 'model' should contain a Fooof-Model by now which will be subsequently saved        #
        ####################################################################################################
        
        #_____ Get FOOOF content _____
        
        # if list -> then multiple fittings := StepWise Fooof
        if isinstance(model,list) == True:
            # get all fields
            fields = ['aperiodic_mode','aperiodic_params_','error_','fooofed_spectrum_','freq_range','freq_res','freqs',
                      'max_n_peaks','min_peak_height','peak_params_','peak_threshold','peak_width_limits','power_spectrum',
                      'r_squared_','gaussian_params_','_peak_fit','_ap_fit','_spectrum_flat','_spectrum_peak_rm']
            #store in add_dict
            add_dict = {}
            # squish all submodels into a single list
            for k in range(len(fields)):
                #collect variable values
                c = []
                for l in range(len(model)):
                    #get model values
                    c.append( eval('model[l].' + fields[k]) )
                #store in dict
                fieldname = fields[k]
                if fieldname[0] == '_': fieldname = fieldname[1:]
                #store data
                add_dict.update( { fieldname: c} )
        
        # if no list -> then single fitting := Fixed / Knee / Initial (bounded offset)
        else:
            #Store content in a dictionary
            add_dict = {'aperiodic_mode': model.aperiodic_mode,
                        'aperiodic_params_': model.aperiodic_params_,
                        'error_': model.error_,
                        'fooofed_spectrum_': model.fooofed_spectrum_,             #fooofed spectrum = fm._bg_fit + fm._peak_fit
                        'freq_range': model.freq_range,
                        'freq_res': model.freq_res,
                        'freqs': model.freqs,
                        'max_n_peaks': model.max_n_peaks,
                        'min_peak_height': model.min_peak_height,
                        'peak_params_': model.peak_params_,
                        'peak_threshold': model.peak_threshold,
                        'peak_width_limits': model.peak_width_limits,
                        'power_spectrum': model.power_spectrum,                  #the original power spectrum (log10-transform)
                        'r_squared_': model.r_squared_,
                        'gaussian_params_': model.gaussian_params_,
                        'peak_fit': model._peak_fit,
                        'ap_fit': model._ap_fit,
                        'spectrum_flat': model._spectrum_flat,                   #spectrum flat = fm.power_spectrum - fm._ap_fit (the aperiodic component)
                        'spectrum_peak_rm': model._spectrum_peak_rm              #spectrum_peak_removed = fm.power_spectrum - fm._peak_fit (the periodic component)
                        }
        
        #add all dict updates in one dict
        all_dict.update({labels[i,0][0]: add_dict})
        #END spectrum loop
        
    #Parcel specific dictionary to save and import in MATLAB
    Parcel_dict = {'Parcel': all_dict}
    
    #Export results to .mat
    sio.savemat(mpath + '/fooof_flexible/' + sub + '/results_fooof/parcel/' + sub + '_parcel_fooof.mat',Parcel_dict)
    #END subject loop
    
#END Script