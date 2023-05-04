%% Script for Correlation Calculation on Parcels (CoM)

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

%load spectra (with properly merged subjects)
load([mpath,'/parcel/ana_power/fooof_spectra.mat']);
%load CoM
load([mpath,'/parcel/ana_CoM/CoM.mat']); CoM = tmp; clear tmp;

%catch problematic subjects
catch_sub = [];

%load cortical grid -> template_grid
load([mpath,'/cortical_grid.mat']);
template_grid.pos = cortical_grid' .* 100; %units: m -> cm
template_grid.unit = 'cm';
clear cortical_grid

%load parcellation
load([mpath,'/parcellation.mat']);

%load significant cluster parcels
load([mpath,'/parcel/ana_stat/Significant_Regions.mat']);

%Parcels that were standing out significantly by the CoE Permutation Test
parcels_of_interest = parcel.masklabel( logical( SigR.hcVSaps .* SigR.pdVSaps ) );

%tests to use [which are in the xlsx-file on seperate sheets]. If two tests are given: First 'test' specifies which test to correlate.
use_test = {'tulia'};%{'updrs_pre_off','tulia','goldenberg_handposition','goldenberg_fingerposition','goldenberg','moca'};
use_field = {'SUMME'};

%updrs      - sum / SUMME
%MoCA       - relative_sum
%goldenberg - SUMME
%tulia      - SUMME

%exclude subjects (py setting CoM.indices to 0 at the specific positions)
if strcmpi(use_test,'moca')
    exclude = {'cbs14'};
    idx = ismember(CoM.sub,exclude);
    gr = fieldnames(CoM.indices);
    for k = 1:length(gr)
        CoM.indices.(gr{k})(idx) = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ Neuropsychological Scores ___%%% -> [more affected bodyside, relative severity, absolute scores]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to xlsx-sheet
testpath = 'C:/Users/mkroe/Documents/Ausbildung/Promotion/project/RestingState/Excel_Subject_Informations/resting_state_study_I/';

CoM.indices.parkinsonism = logical(CoM.indices.cbs);

%groups
groups = {'cbs'};

%excel structure
excel = struct();

for k = 1:length(use_test)
    excel.(use_test{k}) = readtable([testpath,'COHORT_Tests.xlsx'],'Sheet',use_test{k});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ Correlations on Brain | Region of Interest Average___%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%color
Cstat = cbrewer('div','PRGn',64);

%Parcels of Interest Region
idx_ParOIn = ismember(CoM.labels,parcels_of_interest);

%__ I Corr: CoM & Behavioral Scores
p_val_CoE = [];
for gr = 1:length(groups)
    
    %group indices
    gr_indices = CoM.indices.(groups{gr});
    
    %group subjects
    gr_sub = CoM.sub(gr_indices);
    
    %check if test scores exist for subjects
    if sum( ismember(gr_sub,excel.(use_test{1}).subject_id) ) == 0
        warning(['No test scores exist for group ',groups{gr}])
        continue
    end
    
    %get test scores for subjects
    testscores = [];
    for k = 1:length(gr_sub)
        value = excel.(use_test{1}).(use_field{:})(ismember(excel.(use_test{1}).subject_id,gr_sub(k)));
        if iscell(value)
            value = value{:};
        end
        if isstr(value)
            value = str2num(value);
        end
        if isempty(value)
            value = nan;
        end
        testscores = vertcat(testscores,value);
    end
    
    %merge tests (have been normalized to maximally achievable sum before in excel)
    testscores = mean(testscores,2);
    
    if all(isnan(testscores))
        continue
    end
    
    %extract group CoM
    gr_CoM = mean( CoM.CoM( idx_ParOIn, gr_indices ), 1 )';
    
    if any(isnan(testscores))
        gr_CoM = gr_CoM( ~isnan(testscores) );
        testscores = testscores( ~isnan(testscores) );
    end
    
    %Correlation: Calculation
    [Crr,CoM_pval] = corr(gr_CoM,testscores,'Type','Spearman','rows','complete');
    
    p_val_CoE = vertcat(p_val_CoE,CoM_pval);
    
    %bootstrap that stuff for confidence!
    scorr = @(gr_CoM,testscores)(corr(gr_CoM,testscores));
    bootstat = bootstrp(500,scorr,gr_CoM,testscores);
    %histogram of that bootstrap ...  to see if we are actually confident
    figure
    [n_bootstat,x1] = histcounts(bootstat,-1:0.1:1);
    bar(x1(1:end-1) + 0.05,n_bootstat,'histc');
    xlim([-1,1])
    
    %shuffle data to get permutation statistic
    PCrr = [];
    for p = 1:1000
        idx = randperm(length(gr_CoM));
        [PCrr(p),~] = corr(gr_CoM(idx),testscores,'Type','Pearson','rows','complete');
    end
    PCrr = sort(PCrr);
    %Z-Score our empirical correlation value
    Z_Crr = (Crr - mean(PCrr))/std(PCrr);
    %histogram of that permutation distribution ...  how far we are off
    figure
    [n_PCrr,x1] = histcounts(PCrr,-1:0.05:1);
    bar(x1(1:end-1) + 0.025,n_PCrr,'histc');
    hold on
    line([Crr,Crr],[0,max(n_PCrr)],'Color','red','LineWidth',3)
    xlim([-1,1])
    title(['PermStat: CoE | ZScore: ' num2str(Z_Crr), ' Gr.: ',replace(groups{gr},'_',' ')])
    
    figure
    scatter(gr_CoM,testscores,'MarkerFaceColor','k')
    title(['Gr.: ',replace(groups{gr},'_',' '),' | Parcels of Interest | Corr: ',num2str(Crr) ]);
    xlabel('CoE');
    ylabel(['Test Scores']);
    
end
