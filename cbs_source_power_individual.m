%% Script for individual power spectra

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
scp_path = [mpath,'/scripts'];                         %script path

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path);
ft_defaults;

%load cbs project infos, cortex-areal sensor infos, define condition
load([mpath,'/cbs_info.mat']);  %cbs patients info

%load SpecFeat -> only for subjects
load([mpath,'/parcel/ana_shift/spectral_features.mat'])

%load significant cluster parcels
load([mpath,'/parcel/ana_stat/Significant_Regions.mat']);

%load frequency band data
load([mpath,'/parcel/ana_power/fooof_spectra.mat'],'Fooof');

%load parcellation
load([mpath,'/parcellation.mat']);
%load source dataset for creating template grid
load([mpath,'/source/power/hc01/hc01_source_pow.mat']);

%catch problematic subjects
catch_sub = {};

%load cortical grid -> template_grid
load([mpath,'/cortical_grid.mat']);
template_grid.pos = cortical_grid' .* 100; %units: m -> cm
template_grid.unit = 'cm';
%for source investigation
template_grid.labels = psd_for_fooof.label;
%for parcel investigation 
template_grid.masklabel = parcel.masklabel;
template_grid.mask = parcel.mask;
clear cortical_grid psd_for_fooof

%%%%%%%%%%%%%%%%%%
%%%___script___%%%
%%%%%%%%%%%%%%%%%%

%tremor subjects (from pd)
tremor = cbs_clean_subjects(SpecFeat.subjects,cbs_info,{'rest','tremor'},'yes');
%notremor subjects
notremor = cbs_clean_subjects(SpecFeat.subjects,cbs_info,{'rest','tremor'},'-yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ Neuropsychological Scores ___%%% -> [more affected bodyside, relative severity, absolute scores]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to xlsx-sheet
testpath = 'C:/Users/mkroe/Documents/Ausbildung/Promotion/project/RestingState/Excel_Subject_Informations/resting_state_study_I/';

%tests to use [which are in the xlsx-file on seperate sheets]
use_test = {'updrs_pre_off'};%{'updrs_pre_off','tulia','goldenberg_handposition','goldenberg_fingerposition','goldenberg','moca'};
use_field = {'sum'};

%excel structure
excel = struct();

for k = 1:length(use_test)
    excel.(use_test{k}) = readtable([testpath,'COHORT_Tests.xlsx'],'Sheet',use_test{k});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ Set Up Subjects___%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FreqsToUse = ismember(Fooof.freqs, [2:0.5:30]); %Delta: [2 3.5] Theta: [4:7], Alpha: [7.5:13], Beta: [13.5:30],
%subject groups structure
indices.hc = contains(SpecFeat.subjects,'hc');
indices.cbs = contains(SpecFeat.subjects,'cbs');
indices.psp = contains(SpecFeat.subjects,'psp');
indices.aps = logical( contains(SpecFeat.subjects,'cbs') + contains(SpecFeat.subjects,'psp') );
indices.pd = logical( contains(SpecFeat.subjects,'pd') .* ismember(SpecFeat.subjects,notremor));

if strcmpi(use_test,'moca')
    exclude = {'cbs14'};
    idx_exclude = ismember(SpecFeat.subjects,exclude);
    indices.cbs = logical( indices.cbs .* ~idx_exclude );
    indices.aps = logical( indices.aps .* ~idx_exclude );
end

%use overlapping clusters of hc vs. atyp & pd vs. atyp
parcels_of_interest = parcel.masklabel( logical( SigR.hcVSaps .* SigR.pdVSaps ) );

%___ Parcel Spectra Plots

%groupnames
groupnames = fieldnames(indices);

%code group comparison [row vs column]
disp( strcat(num2str((1:length(groupnames))'),': ',groupnames) )
%use indices to refer to groups
CodeA = [1 4 5];
%group colors
col = cbrewer('qual', 'Set3', 12, 'pchip');
col = col([1,4,6],:);

figure('Position',[10 10 1500 450])
for n = 1:length(CodeA)
    subplot(1,3,n)
    %group spectra
    code = CodeA(n);
    %ROI indices
    code_ROI = find( contains(Fooof.labels,parcels_of_interest) );
    %use corresponding indices
    specs = squeeze( Fooof.fooof_spec(code_ROI,:,indices.(groupnames{code})) );
    %group average
    gr_avg = squeeze( mean(specs,1) );
    %plot spectra
    hold on
    for s = 1:size(specs,1)
        plot( Fooof.freqs, squeeze(specs(s,:,:)),'Color',col(n,:),'LineWidth', .1)
    end
    %plot( Fooof.freqs, gr_avg,'Color','k','LineWidth',0.5);
    %title(replace(upper(groupnames{CodeA(n)}),'_','. '))
    xlim([Fooof.freqs(1),Fooof.freqs(end)])
    ylim([-0.1,1.1])
    xlabel('Frequency')
    ylabel('log10(Power) [without 1/f]')
    set(gca,'FontSize',17)
    legend({replace(upper(groupnames{CodeA(n)}),'_','. ')},'FontSize',17)
    legend boxoff
end
print([mpath,'/parcel/ana_power/Individual_Parcel_Spectra_in_ROI_with_ROI_avg.tiff'],'-dtiff','-r300');

close all
