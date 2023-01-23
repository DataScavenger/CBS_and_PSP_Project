%% cbs parcel connectivity analysis [source level] & parcel power analysis

%clean workspace and command window
clearvars;
clc;

%%%%%%%%%%%%%%%%%%%%
%%%___settings___%%%
%%%%%%%%%%%%%%%%%%%%

%path settings
mpath = 'C:/data';                                     %mainpath
ft_path = 'C:/toolboxes/fieldtrip-20201214';           %fieltrip path
fct_path = [mpath,'/functions'];                       %function path (my own functions)
roi_path = 'C:/toolboxes/MEG-ROI-nets-master';         %orthogonalization
scp_path = [mpath,'/scripts'];                         %script path

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path,roi_path);
ft_defaults;

%load cbs info
load([mpath,'/cbs_info.mat']);  %cbs patients info
subjects = fieldnames(cbs_info);

%load parcellation
load([mpath,'/parcellation.mat'])

%catch problematic subjects
catch_sub = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ Parcel Power ___%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(subjects)
    
    try
        %load parcel time-courses
        load([mpath,'/parcel/time/',subjects{i},'/',subjects{i},'_parcel_time.mat'])
        
        if ~exist([mpath,'/parcel/power/',subjects{i}],'dir')
            mkdir([mpath,'/parcel/power/',subjects{i}]);
        end
        
        %transfer to frequency domain & compute csd
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'fourier';
        cfg.channel = 'all';
        cfg.taper = 'dpss';
        cfg.foilim = [1 48];
        cfg.tapsmofrq = 2;
        parcel_freq = ft_freqanalysis(cfg,parcel_time);
        
        %get power
        cfg = [];
        psd_for_fooof = ft_freqdescriptives(cfg,parcel_freq);
        %clean the structure
        psd_for_fooof = rmfield(psd_for_fooof,{'dimord','cumtapcnt','cumsumcnt','cfg'});
        
        %save psd
        save([mpath,'/parcel/power/',subjects{i},'/',subjects{i},'_parcel_pow.mat'],'psd_for_fooof');
    catch
        catch_sub = horzcat(catch_sub,i);
    end
    
end
