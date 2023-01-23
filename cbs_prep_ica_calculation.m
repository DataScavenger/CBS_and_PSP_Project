%% Independent Component Analysis

%clean workspace and command window
clearvars;
clc;

%path settings
mpath = 'C:/data/';                                    %mainpath
ft_path = 'C:/toolboxes/fieldtrip-20201214';           %fieltrip path
fct_path = [mpath,'/functions'];                       %function path (my own functions)
scp_path = [mpath,'/scripts'];                         %script pathft_defaults;

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path);
ft_defaults;

%load cbs project infos, cortex-areal sensor infos, define condition
load([mpath,'/cbs_info.mat']);  %cbs patients info
subjects = fieldnames(cbs_info);
condition = 'rest';

%Calcuate the ICA on gradiometers
use_sensors = 'seperately';
these = []; % datasets did not work. Shows index i
v = 1;

disp( strcat(num2str((1:length(subjects))'),' :',subjects) )
i = find(contains(subjects,'pd65')); i=i(1);

for i = 1:length(subjects)
    
    try
    
    components = []; %Structure to put in independent components
    
    %% prelimnaries: preprocessing | ICA
    
    %filepath
    filepath = cbs_info.(subjects{i}).(condition).path;
    if iscell(filepath);    filepath = filepath{1};     end
    
    %load data
    cfg = [];
    cfg.dataset = [filepath,'.fif'];
    data = ft_preprocessing(cfg);   
    
    %delete bad channels
    cfg = [];
    cfg.channel =  ['meg',cbs_info.(subjects{i}).(condition).badchan(:)'];
    data_select = ft_selectdata(cfg,data);
    
    %downsample before hp filter to 1000Hz (if necessary) --> High Sampling Rate as recommended by Steve Luck (See Solutiontoeverything.txt)
    if data_select.fsample > 1000; cfg = []; cfg.resamplefs = 1000; cfg.detrend = 'no'; data_select = ft_resampledata(cfg,data_select); end
    
    %hp filter
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfiltertype = 'fir';
    cfg.hpfreq = 1;
    data_select = ft_preprocessing(cfg,data_select);
    
    %lp-filter
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfiltertype = 'fir';
    cfg.lpfreq = 42;
    data_select = ft_preprocessing(cfg,data_select);
    
    %dummy trials for ica [to delete segments]
    cfg = [];
    cfg.length = 1;
    cfg.overlap = 0;
    data_select_seg = ft_redefinetrial(cfg,data_select);
    
    %load the artifactual sample information    
    load(['C:/data/artifacts/',subjects{i},'_rest.mat'])
    
    timepoints_to_delete = [];
    %translate artifact structure data from 250Hz resampled data (for artifacts) to 1000Hz resampled data
    if isfield(ArtiMEG.artfctdef,'visual') && isfield(ArtiMEG.artfctdef.visual,'artifact')
        ArtiMEG.artfctdef.visual.artifact = ArtiMEG.artfctdef.visual.artifact .* 4;
        timepoints_to_delete = [ timepoints_to_delete ; ArtiMEG.artfctdef.visual.artifact ]; 
    end
    if isfield(ArtiMEG.artfctdef,'spectra') && isfield(ArtiMEG.artfctdef.spectra,'artifact')
        ArtiMEG.artfctdef.spectra.artifact = ArtiMEG.artfctdef.spectra.artifact .* 4;
        timepoints_to_delete = [ timepoints_to_delete ; ArtiMEG.artfctdef.spectra.artifact ]; 
    end
    if isfield(ArtiMEG.artfctdef,'emg') && isfield(ArtiMEG.artfctdef.emg,'artifact')
        ArtiMEG.artfctdef.emg.artifact = ArtiMEG.artfctdef.emg.artifact .* 4;
        timepoints_to_delete = [ timepoints_to_delete ; ArtiMEG.artfctdef.emg.artifact ];
    end
    if isfield(ArtiMEG.artfctdef,'protocoll') && isfield(ArtiMEG.artfctdef.protocoll,'artifact')
        ArtiMEG.artfctdef.protocoll.artifact = ArtiMEG.artfctdef.protocoll.artifact .* 4;
        timepoints_to_delete = [ timepoints_to_delete ; ArtiMEG.artfctdef.protocoll.artifact ];
    end
    %Note: I did not include the field: .tremor -> So a person showing both periods of tremor and no tremor, or things in between this kind of tremor mixed data is forwarded to the ICA. Distinguish tremor versus no-tremor periods with this field in scripts later.
    
    %BE SURE THAT ARTIFACTS WERE ACTUALLY LOOKED AT WITH 250HZ
    
    %Delete artifacts
    cfg = [];
    cfg.artfctdef.visual.artifact = timepoints_to_delete;
    data_clean = ft_rejectartifact(cfg,data_select_seg);
    
    %get matrix rank
    if rank(cell2mat(data_clean.trial)) ~= length(data_clean.label)  
       error('Deficiency of Matrix Rank') 
    end
    
    %% ICA on gradiometers
    
    if strcmpi(use_sensors,'seperately')
        %make ICA for gradiometers
        cfg = [];
        cfg.channel = {'megplanar'};
        gradio_data = ft_selectdata(cfg,data_clean);
        
        %run ica
        cfg = [];
        cfg.method = 'runica';
        gradio_comp = ft_componentanalysis(cfg, gradio_data);
        
        %save ica results
        components.gradio_components = gradio_comp;
        
        %Save ica results
        if ~exist(['D:\more_clean_data\ica\',subjects{i}],'dir')
            mkdir(['D:\more_clean_data\ica\',subjects{i}])
        end
        save(['D:\more_clean_data\ica\',subjects{i},'\',subjects{i},'_independent_components.mat'],'components');
        
    end
                
    catch
        these(v) = i;
        v = v+1;
    end
end