%% Visualize the calculated independent components and note down the bad components.

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
load([mpath,'/cbs_info.mat']);  %subjects/patients info
subjects = fieldnames(cbs_info);
condition = {'rest'};
which_ICs = {'independent_components'};

%the subject number to enter in load( ... subject{choice} ... )
disp( strcat(num2str((1:length(subjects))'),' :',subjects) )
choice = 135;
disp(subjects{choice});

%load components of subject -> variable name 'components'
load(['D:/more_clean_data/ica/',subjects{choice},'/',subjects{choice},'_independent_components.mat']);

%choose gradiometers / magnetometers for independent component selection
sensor_type = fieldnames(components);
disp( strcat(num2str((1:length(sensor_type))'),' :',sensor_type) );
choice_sensor_type = 1; %this chooses gradiometers
%use these components
comp = components.(sensor_type{choice_sensor_type});

%visualize components
cfg = [];
cfg.channel = 1:5;
if contains(sensor_type{choice_sensor_type},'gradio'); cfg.layout = 'neuromag306planar.lay'; end
if contains(sensor_type{choice_sensor_type},'magneto'); cfg.layout = 'neuromag306mag.lay'; end
cfg.compscale  = 'local';
cfg.continuous = 'yes';
if contains(sensor_type{choice_sensor_type},'gradio'); cfg.ylim = [-1e-11,1e-11]; end
if contains(sensor_type{choice_sensor_type},'magneto'); cfg.ylim = [-1e-14,1e-14]; end
ft_databrowser(cfg, comp);
