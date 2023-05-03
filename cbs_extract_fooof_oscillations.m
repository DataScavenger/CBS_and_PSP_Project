%% function to extract oscillatory peaks from fooof flexible ['fixed','knee','initial','stepwise']

% FooofData:    Output of cbs_parcel_fooof_ana_single_flexible.py

function [ oscillations,amp ] = cbs_extract_fooof_oscillations(FooofData)
    
    oscillations = [];
    amp = [];
    
    % depending on the fitting - modus different code must be executed
    switch FooofData.aperiodic_mode(1,:)
        
        case 'fixed'
            
            if ~isempty( FooofData.gaussian_params_ )
                % get oscillations
                oscillations = FooofData.gaussian_params_(:,1);
                % get amplitude
                amp = FooofData.gaussian_params_(:,2);
            end
                        
        case 'knee'
            
            if ~isempty( FooofData.gaussian_params_ )
                % get oscillations
                oscillations = FooofData.gaussian_params_(:,1);
                % get amplitude
                amp = FooofData.gaussian_params_(:,2);
            end
            
        case 'initial'
            
            if ~isempty( FooofData.gaussian_params_ )
                % freq resolution
                freqres = FooofData.freqs(2) - FooofData.freqs(1);
                % get oscillations
                oscillations = ( freqres * FooofData.gaussian_params_(:,1) ) + FooofData.freqs(1) - freqres;
                % get amplitude
                amp = FooofData.gaussian_params_(:,2);
            end
            
        case 'stepwise'
            
            for k = 1:length(FooofData.freqs)
                if iscell(FooofData.gaussian_params_)
                    
                    if ~isempty( FooofData.gaussian_params_{k} )
                        % freq resolution
                        freqres = FooofData.freqs{k}(2) - FooofData.freqs{k}(1);
                        % get oscillations
                        oscillations =  vertcat( oscillations, ( freqres * FooofData.gaussian_params_{k}(:,1) ) + FooofData.freqs{k}(1) - freqres );
                        % get amplitude
                        amp = vertcat( amp, FooofData.gaussian_params_{k}(:,2) );
                    end
                    
                else
                    %this is because in very rare cases python does not put the gaussian params into a cell, but into a matrix
                    if ~isempty( FooofData.gaussian_params_ )
                        % freq resolution
                        freqres = FooofData.freqs{k}(2) - FooofData.freqs{k}(1);
                        % get oscillations
                        oscillations =  vertcat( oscillations, ( freqres * FooofData.gaussian_params_(k,:,1) ) + FooofData.freqs{k}(1) - freqres );
                        % get amplitude
                        amp = vertcat( amp, FooofData.gaussian_params_(k,:,2) );
                    end
                    
                end
                
            end
    end
    
end