%% function to extract from power fooof flexible ['fixed','knee','initial','stepwise']

% FooofData:        Output of cbs_parcel_fooof_ana_single_flexible.py
% choose_spectrum:  fooofed_spectrum or power_spectrum
% freq_range:       Frequency range to fit the spectral data in 'FooofData' to
% freq_res:         Frequency resolution of the FooofData

function [spectrum,freqs] = cbs_extract_fooof_flexible(FooofData,choose_spectrum,freq_range,freq_res)
    
    % some parameters [freqbins and spectrum]
    freqbins    = freq_range(1):freq_res:freq_range(2);
    spectrum    = nan(1,length(freqbins));
    
    % ensure that freqbins contains als frequencies that must be covered as given by FooofData
    if iscell(FooofData.freqs)
       FreqCheck = [ FooofData.freqs{:} ];
    else
       FreqCheck = FooofData.freqs;
    end
    
    if ~( all( ismember( FreqCheck, freqbins ) ) )
        error('"freq_range" or "freq_res" are not fitting with the information from "FooofData".')
    end
    
    % depending on the fitting - modus different code must be executed
    switch FooofData.aperiodic_mode(1,:)
        
        case 'fixed'
            % get fooof spectrum
            fooof_spectrum = FooofData.(choose_spectrum);
            % prune fooof spectrum
            spectrum( ismember( freqbins, FooofData.freqs ) ) = fooof_spectrum;
            
            freqs = sum(freqbins,1,'omitnan');
            
        case 'knee'
            % get fooof spectrum
            fooof_spectrum = FooofData.(choose_spectrum);
            % prune fooof spectrum
            spectrum( ismember( freqbins, FooofData.freqs ) ) = fooof_spectrum;
            
            freqs = sum(freqbins,1,'omitnan');
            
        case 'initial'
            % get fooof spectrum
            fooof_spectrum = FooofData.(choose_spectrum);
            % prune fooof spectrum
            spectrum( ismember( freqbins, FooofData.freqs ) ) = fooof_spectrum;
            
            freqs = sum(freqbins,1,'omitnan');
            
        case 'stepwise'
            % get fooof spectra
            fooof_freqs = nan( length(FooofData.freqs), length(freqbins) );
            fooof_spectra = nan( length(FooofData.freqs), length(freqbins) );
            
            for k = 1:length(FooofData.freqs)
                fooof_freqs(k, ismember( freqbins, FooofData.freqs{k} ) )   = FooofData.freqs{k};
                fooof_spectra(k, ismember( freqbins, FooofData.freqs{k} ) ) = FooofData.(choose_spectrum){k};
            end
            
            overlap = find( ~isnan( sum(fooof_freqs,1) ) );
            
            for k = 1:length(overlap)
                fooof_freqs(k,overlap) = nan;
                fooof_spectra(k,overlap) = nan;
            end
            
            freqs = sum(fooof_freqs,1,'omitnan'); freqs(freqs == 0) = nan;
            spectrum = sum(fooof_spectra,1,'omitnan'); spectrum(isnan(freqs)) = nan;
            
    end
    
end