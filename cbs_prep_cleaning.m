%% Clean Data Script

%clean workspace and command window
clearvars;
clc;

%% path settings
mpath = 'C:/data';                                     %mainpath
ft_path = 'C:/toolboxes/fieldtrip-20201214';           %fieltrip path
fct_path = [mpath,'/functions'];                       %function path (my own functions)
scp_path = [mpath,'/scripts'];                         %script path
cbs_info_path = [mpath,'/cbs_info.mat'];               %cbs_info path
csinfo_path = [mpath,'/csinfo.mat'];                   %csinfo path

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path);
ft_defaults;

%load cbs project infos, cortex-areal sensor infos, define condition
load(cbs_info_path); %cbs patients info
subjects = fieldnames(cbs_info);

load(csinfo_path); %for power spectra inspection

%create folder for clean data
if ~exist('D:/more_clean_data/clean_data/','dir')
    mkdir('D:/more_clean_data/clean_data/');
end

%% script parameters
condition = 'rest';
%delete or keep presumably bad trials ('yes','no')
delete_bad = 'yes';

%Use independent components ('gradio_components')
use_ic = {'gradio_components'};
%Load indices of bad components
load('D:/more_clean_data/ica/bad_independent_components.mat')

%These are the subjects
disp( strcat(num2str((1:length(subjects))'),' :',subjects) )
%choose subject
i = find(strcmpi(subjects,'pd65'));

%index to look for datasets with error
catch_subjects = [];

processed_subjects =  [];

for i = 1:length(subjects)
    
    try
        %% Time-Domain Clean Data
        
        %raw data: filepath
        filepath = cbs_info.(subjects{i}).(condition).path;
        if iscell(filepath);    filepath = filepath{1};     end
        
        %load data
        cfg = [];
        cfg.dataset = [filepath,'.fif'];
        data = ft_preprocessing(cfg);
        
        %delete bad channels
        cfg = [];
        cfg.channel = ['meg',cbs_info.(subjects{i}).(condition).badchan(:)'];
        data_select = ft_selectdata(cfg,data);
        
        %downsample before hp filter to 1000Hz
        if data_select.fsample > 1000; cfg = []; cfg.resamplefs = 1000; cfg.detrend = 'no'; data_select = ft_resampledata(cfg,data_select); end
        
        %hp filter
        cfg = [];
        cfg.hpfilter = 'yes';
        cfg.hpfiltertype = 'fir';
        cfg.hpfreq = 1;
        data_select = ft_preprocessing(cfg,data_select);
        
        if ~isempty(use_ic) && strcmpi(use_ic,'gradio_components')
            
            %load independent components            
            load(['D:/more_clean_data/ica/',subjects{i},'/',subjects{i},'_independent_components.mat'])
            
            %data only gradiometers
            cfg = [];
            cfg.channel = 'megplanar';
            data_select_gradios = ft_selectdata(cfg,data_select);
            
            %reject gradiometer components from gradiometers
            cfg = []; %cfg.demean = 'yes' per default
            cfg.updatesens = 'yes';
            cfg.component = bad_independent_components.(subjects{i}).(condition).gradio_components; %select indices of bad components
            data_select_clean = ft_rejectcomponent(cfg, components.gradio_components, data_select_gradios);
            
        else
            
            data_select_clean = data_select; %data_select_clean
            
        end
        
        %downsample
        cfg = [];
        cfg.resamplefs = 250;
        cfg.detrend = 'no';
        data_select_clean = ft_resampledata(cfg, data_select_clean);
        
        %load artifactual information
        load(['D:/more_clean_data/artifact/',subjects{i},'_rest.mat'])
        
        %If the artifactual sample information exists, then use it for cleaning
        if exist('ArtiMEG','var')
            % collect bad trials from artifact fields (information that has been collected in cbs_prep_databrowse.m)
            timepoints_to_delete = []; artfields = {'visual','spectra','emg','eyesclosed','tremor','protocoll'};
            for k = 1:length(artfields); if isfield(ArtiMEG.artfctdef,artfields{k}); timepoints_to_delete = [ timepoints_to_delete; ArtiMEG.artfctdef.(artfields{k}).artifact ]; end; end
            
            cfg = [];
            cfg.artfctdef.reject = 'partial';
            cfg.artfctdef.visual.artifact = timepoints_to_delete;
            data_select_clean = ft_rejectartifact(cfg,data_select_clean);
        end
        
        %save the data
        save(['D:/more_clean_data/clean_data/',subjects{i},'_',condition,'_time.mat'],'data_select_clean');
        
    catch
        
        catch_subjects = horzcat(catch_subjects,i);
        
    end
    
end