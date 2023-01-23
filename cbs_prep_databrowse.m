%% Browse data to find bad channels & bad time segments

%clean workspace and command window
clearvars;
clc;

%path settings
mpath = 'C:/data';                                     %mainpath
ft_path = 'C:/toolboxes/fieldtrip-20201214';           %fieltrip path
fct_path = [mpath,'/functions'];                       %function path (my own functions)
scp_path = [mpath,'/scripts'];                         %script path

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path);
ft_defaults;

%load cbs project infos, cortex-areal sensor infos, define condition
load('C:/data/cbs_info.mat');              %structure -> with field 'subject'-ID -> with subfields filepath, birthday, ...
load('C:/data/csinfo.mat');                %for power spectra inspection
subjects = fieldnames(cbs_info);

condition = 'rest';
check_spectra = 'yes';

                      %structure    . field       . subfield
%load protocoll times (protocoll -> subject_id -> timepoints (sec.) to delete given by experimental protocoll)
load([mpath,'/artifacts/protocoll_times.mat'])

%These are the subjects
disp( strcat(num2str((1:length(subjects))'),' :',subjects) )
%Run this if you want to choose a specific subject 'e.g.pd22'
i = find(contains(subjects,'pd65')); i=i(1);

for i = 1:length(subjects)
    
    %filepath to raw data
    filepath = cbs_info.(subjects{i}).(condition).path;
    if iscell(filepath);    filepath = filepath{1};     end
    %Load raw data
    cfg = [];
    cfg.dataset = [filepath,'.fif'];
    data = ft_preprocessing(cfg);
    
    %downsample
    cfg = [];
    cfg.resamplefs = 250;
    cfg.detrend = 'no';
    data = ft_resampledata(cfg, data);
    
    %hp-filter
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfiltertype = 'fir';
    cfg.hpfreq = 0.5;
    data = ft_preprocessing(cfg,data);
    
    %load subject protocoll (in seconds -> samples)
    if isfield(Protocoll,subjects{i}) && ~isempty( Protocoll.(subjects{i}).protocoll )
        subproto_sec = Protocoll.(subjects{i}).protocoll;
        subproto = subproto_sec * data.fsample;
        if subproto(1) == 0; subproto(1) = 1; end
    else
        subproto = double.empty(0,2);
        subproto_sec = [];
        warning(['The subject ',subjects{i},' is not within the Protocoll-structure'])
    end
    
    %look at data manually [first EOG only]
    cfg = [];
    cfg.channel = {'EOG*'};
    cfg.colorgroups = 'allblack';
    cfg.viewmode = 'vertical';
    cfg.ylim = [-0.0002,0.0002];
    ArtiEOG = ft_databrowser(cfg,data);
    
    %run this after EOG inspection
    blinks = ArtiEOG.artfctdef.visual.artifact;     %Will be disregarded for ICA anyways
    eyesclosed = ArtiEOG.artfctdef.visual.artifact; %Suspicious phases of absent eye activity? E.g. if patient reported to be tired
    
    %% magnetometer inspection
    cfg = [];
    cfg.channel = 'MEG***1';
    tmp = ft_selectdata(cfg,data);
    cfg = [];
    %optional fields
    cfg.artfctdef.blinks.artifact = blinks;
    cfg.artfctdef.eyesclosed.artifact = eyesclosed;
    cfg.artfctdef.protocoll.artifact = subproto;
    cfg.artfctdef.visual.artifact = double.empty(0,2);
    cfg.channel = {'MEG***1'};
    cfg.colorgroups = 'allblack';
    cfg.viewmode = 'vertical';
    ArtiMEG = ft_databrowser(cfg,tmp);
    
    %compute power spectrum (magnetometers)
    if strcmp(check_spectra,'yes')
        
        %segmented data 1 sec
        cfg = [];
        cfg.length = 1;
        cfg.overlap = 0;
        tmp_segments = ft_redefinetrial(cfg,tmp);
        
        %Fourier transformation and extraction of power
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.output = 'pow';
        cfg.foilim = [1 40];
        cfg.keeptrials = 'yes';
        tmp_segments_pow = ft_freqanalysis(cfg,tmp_segments);
                
        %area names                    %indexing only the fields with regional information
        cs_names = fieldnames(csinfo); num_areas = find(strncmpi(cs_names,'MEG',3)); num_areas = num_areas(1) - 1;
        cs_names = cs_names(1:num_areas);
        %corresponding data
        area_info = struct2cell(csinfo);
        area_info = area_info(1:num_areas);
        
        %clean area_info from gradiometers
        for x = 1:length(area_info);    area_info{x} = area_info{x}( endsWith(area_info{x},'1') );    end
        
        %average data & plot power spectra
        for x = 1:length(cs_names)
            
            %average: area data to append
            cfg = [];
            cfg.channel = area_info{x};
            cfg.avgoverchan = 'yes';
            tmp_segments_areapow = ft_selectdata(cfg,tmp_segments_pow);
            %change to log10 normalization
            tmp_segments_areapow.powspctrm = log10( tmp_segments_areapow.powspctrm );
            
            %make waterfall
            cfg = [];
            cfg.xlim = [1 40];      %frequency
            cfg.ylim = [1 size(tmp_segments_areapow.powspctrm,1)];       %trials
            cfg.zlabel = 'Power [log10 normalized]';
            figure;
            cbs_waterfall(cfg,squeeze(tmp_segments_areapow.powspctrm));
            region = cs_names{x}; to_be_spaced = strfind(region,'_'); region(to_be_spaced) = ' ';
            title({['subject: ',subjects{i},' | sensor average: ',region]});
        end
        
    end
    %enter spectra indices that should be cleaned
    ArtiMEG.artfctdef.spectra.artifact = cbs_trials2samples([],tmp_segments.sampleinfo);
    
    %% gradiometer inspection
    cfg = [];
    cfg.channel = 'megplanar';
    tmp = ft_selectdata(cfg,data);
    cfg = [];
    cfg.artfctdef = ArtiMEG.artfctdef;
    cfg.colorgroups = 'allblack';
    cfg.viewmode = 'vertical';
    ArtiMEG = ft_databrowser(cfg,tmp);
    
    %timepoints to trials
    timepoints_to_trials = vertcat( ArtiMEG.artfctdef.eyesclosed.artifact / 250, ArtiMEG.artfctdef.protocoll.artifact / 250);
    timepoints_to_trials(:,1) = ceil(timepoints_to_trials(:,1));
    timepoints_to_trials(:,2) = ceil(timepoints_to_trials(:,2));
    badtrials = [];
    for k= 1:size(timepoints_to_trials,1)
        badtrials = [badtrials,timepoints_to_trials(k,1):timepoints_to_trials(k,2)];
    end
    
    %compute power spectrum (gradiometers)
    if strcmp(check_spectra,'yes')
        
        %segmented data 1 sec
        cfg = [];
        cfg.length = 1;
        cfg.overlap = 0;
        tmp_segments = ft_redefinetrial(cfg,tmp);
                
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.output = 'pow';
        cfg.foilim = [1 40];
        cfg.keeptrials = 'yes';
        tmp_segments_pow = ft_freqanalysis(cfg,tmp_segments);
                
        if exist('badtrials','var')
            cfg  = [];
            cfg.trials = 1:size(tmp_segments_pow.powspctrm,1);
            cfg.trials(badtrials) = [];
            tmp_segments_pow = ft_selectdata(cfg,tmp_segments_pow);
        end
        
        %area names                    %indexing only the fields with regional information
        cs_names = fieldnames(csinfo); num_areas = find(strncmpi(cs_names,'MEG',3)); num_areas = num_areas(1) - 1;
        cs_names = cs_names(1:num_areas);
        %corresponding data
        area_info = struct2cell(csinfo);
        area_info = area_info(1:num_areas);
        
        %clean area_info from gradiometers
        for x = 1:length(area_info);    area_info{x}( endsWith(area_info{x},'1') ) = [];    end
        
        %average data & plot power spectra
        for x = 1:length(cs_names)
            
            %average: area data to append
            cfg = [];
            cfg.channel = area_info{x};
            cfg.avgoverchan = 'yes';
            tmp_segments_areapow = ft_selectdata(cfg,tmp_segments_pow);
            %change to log10 normalization
            tmp_segments_areapow.powspctrm = log10( tmp_segments_areapow.powspctrm );
            
            %make waterfall
            cfg = [];
            cfg.xlim = [1 40];         %frequency
            cfg.XTickLabel = {'0','10','20','30','40'};
            cfg.ylim = [1 size(tmp_segments_areapow.powspctrm,1)];       %trials
            cfg.zlabel = 'Power [log10 normalized]';
            figure;
            cbs_waterfall(cfg,squeeze(tmp_segments_areapow.powspctrm));
            region = cs_names{x}; to_be_spaced = strfind(region,'_'); region(to_be_spaced) = ' ';
            title({['subject: ',subjects{i},' | sensor average: ',region]});
            
        end
        
    end
    %enter spectra indices that should be cleaned
    ArtiMEG.artfctdef.spectra.artifact = [ ArtiMEG.artfctdef.spectra.artifact; cbs_trials2samples([],tmp_segments.sampleinfo) ];
    
    %% EMG inspection (esp. for tremor / restlessness in patient datasets)   
    cfg = [];
    cfg.channel = 'EMG*';
    tmp = ft_selectdata(cfg,data);
    
    %segmented data 1 sec and make waterfall
    cfg = [];
    cfg.length = 1;
    cfg.overlap = 0;
    tmp_segments = ft_redefinetrial(cfg,tmp);
    
    %Fourier transformation and power extraction
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.output = 'pow';
    cfg.foilim = [1 40];
    cfg.keeptrials = 'yes';
    tmp_segments_pow = ft_freqanalysis(cfg,tmp_segments);
    
    %create left and right emg power spectra
    cfg = [];
    cfg.channel = {'*62'};
    cfg.avgoverchan = 'yes';
    tmp_segments_leftemg = ft_selectdata(cfg,tmp_segments_pow);
    tmp_segments_leftemg.powspctrm = log10( tmp_segments_leftemg.powspctrm );
    cfg = [];
    cfg.channel = {'*61'};
    cfg.avgoverchan = 'yes';
    tmp_segments_righemg = ft_selectdata(cfg,tmp_segments_pow);
    tmp_segments_righemg.powspctrm = log10( tmp_segments_righemg.powspctrm );
    
    %%% Skip this block if you don't want to delete time segments %%%
    
    %timepoints to trials
    timepoints_to_trials = ArtiMEG.artfctdef.spectra.artifact / 250;
    timepoints_to_trials(:,1) = ceil(timepoints_to_trials(:,1));
    timepoints_to_trials(:,2) = ceil(timepoints_to_trials(:,2));
    %timepoints that have already been marked as bad
    badtrials = [];
    for k= 1:size(timepoints_to_trials,1)
        badtrials = [badtrials,timepoints_to_trials(k,1):timepoints_to_trials(k,2)];
    end
    tmp_segments_leftemg.powspctrm(badtrials,:,:) = nan;
    tmp_segments_righemg.powspctrm(badtrials,:,:) = nan;
    %eventuelly include information from protocoll to delete
    bad_trials = []; for k = 1:size(subproto_sec,1); bad_trials = [ bad_trials, subproto_sec(k,1):subproto_sec(k,2) ]; end
    tmp_segments_leftemg.powspctrm(bad_trials,:,:) = nan;
    tmp_segments_righemg.powspctrm(bad_trials,:,:) = nan;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg = [];
    cfg.xlim = [1 40];         %frequency
    cfg.ylim = [1 size(tmp_segments_leftemg.powspctrm,1)];       %trials
    cfg.zlabel = 'Power [log10 normalized]';
    figure;
    cbs_waterfall(cfg,squeeze(tmp_segments_leftemg.powspctrm));
    title('left')
    cfg = [];
    cfg.xlim = [1 40];         %frequency
    cfg.ylim = [1 size(tmp_segments_righemg.powspctrm,1)];       %trials
    cfg.zlabel = 'Power [log10 normalized]';
    figure;
    cbs_waterfall(cfg,squeeze(tmp_segments_righemg.powspctrm));
    title('right')
    
    %optional fields(include timestamps that should be deleted based on emg power spectra)
    ArtiMEG.artfctdef.tremor.artifact = cbs_trials2samples([],tmp_segments.sampleinfo);
    ArtiMEG.artfctdef.emg.artifact = cbs_trials2samples([],tmp_segments.sampleinfo);
    
    %timecourse emg (run for noisy emg periods OR tremor. Eventually run the code twice)
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfiltertype = 'fir';
    cfg.hpfreq = 10;
    cfg.rectify = 'yes';
    tmp = ft_preprocessing(cfg,tmp);
    
    cfg = [];
    cfg.colorgroups = 'allblack';
    cfg.viewmode = 'vertical';
    cfg.ylim = [-2*10^(-5),2*10^(-5)];
    ArtiMEG_EMG = ft_databrowser(cfg,tmp);
    
    %optional fields (delete timestamps that should be deleted based on emg time series)
    ArtiMEG.artfctdef.emg.artifact = [ ArtiMEG.artfctdef.emg.artifact; ArtiMEG_EMG.artfctdef.visual.artifact ];       %noisy emg periods
    ArtiMEG.artfctdef.tremor.artifact = [ ArtiMEG.artfctdef.tremor.artifact; ArtiMEG_EMG.artfctdef.visual.artifact ]; %tremor periods
    
    %save the artifact structure
    sub = '';
    save(['D:\more_clean_data\artifact\',sub,'_rest.mat'],'ArtiMEG');
end
