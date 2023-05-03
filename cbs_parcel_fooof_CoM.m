%% Script for CoM Calculation and Investigation with Parcels

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
rain_path= 'C:/toolboxes/RainCloudPlots-master/tutorial_matlab';       %raincloud path

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path,rain_path);
ft_defaults;

%load cbs project infos, cortex-areal sensor infos, define condition
load([mpath,'/cbs_info.mat']);  %cbs patients info

%load spectra
load([mpath,'/parcel/ana_power/fooof_spectra.mat']);
%load source information (for positions, unit, labels)
load([mpath,'/parcel/time/hc01/hc01_parcel_time.mat']);

%catch problematic subjects
catch_sub = [];

%load cortical grid -> template_grid
load([mpath,'/cortical_grid.mat']);
template_grid.pos = cortical_grid' .* 100; %units: m -> cm
template_grid.unit = 'cm';
clear cortical_grid

%%%%%%%%%%%%%%%%%%
%%%___script___%%%
%%%%%%%%%%%%%%%%%%

%tremor subjects (from pd)
tremor = cbs_clean_subjects(Fooof.sub,cbs_info,{'rest','tremor'},'yes');
%notremor subjects
notremor = cbs_clean_subjects(Fooof.sub,cbs_info,{'rest','tremor'},'-yes');

%%%%%%%%%%%%%%%%%
%%%___ CoM ___%%%
%%%%%%%%%%%%%%%%%

%Center of Mass ... or Frequency
CoM = zeros( size(Fooof.fooof_spec,1),size(Fooof.fooof_spec,3) );
FreqsToUse = ismember(Fooof.freqs, [4:0.5:30]); %Theta: [4:7], Alpha: [7:13], Beta: [13:30], Mu: [8:20]
%subject groups structure
indices.hc = contains(Fooof.sub,'hc');
indices.cbs = contains(Fooof.sub,'cbs');
indices.psp = contains(Fooof.sub,'psp');
indices.aps = logical( contains(Fooof.sub,'cbs') + contains(Fooof.sub,'psp') );
indices.pd = logical( contains(Fooof.sub,'pd') .* ismember(Fooof.sub,notremor));

%groupnames
groupnames = fieldnames(indices);

%loop parcels and subjects
for p = 1:size(Fooof.fooof_spec,1)
    for s = 1:size(Fooof.sub,1)
        
        %extract weights (log10(power) in this frequency that has been corrected for 1/f estimate (by substraction) -> thus are positive numbers. Output from fooof-algorithm is log10(power) values)
        vec = Fooof.fooof_spec(p,FreqsToUse,s);
        
        %reference to the smallest value in there (will be close to 0)
        vec = vec - min(vec);
        
        %Calculate CoM (sum over weighted frequency bins. Please note the "weights" in vec here are log10(power) values. The weights are all positive, as they have been corrected for the 1/f estimate)
        CoM(p,s) = (vec * Fooof.freqs(FreqsToUse)') / sum(vec);
        
    end
end

%variable settings
tmp.CoM = CoM;
tmp.labels = parcel_time.label;
tmp.sub = Fooof.sub;
tmp.indices = indices;

%path to save figures
if ~exist([mpath,'/parcel/ana_CoM/'],'dir')
    mkdir([mpath,'/parcel/ana_CoM/']);
end
save([mpath,'/parcel/ana_CoM/CoM.mat'],'tmp');

clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ CoM T-Values & Plots ___%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mesh in mm
mesh = ft_read_headshape([ft_path,'/template/anatomy/surface_white_both.mat']);

%put t-values on the brain
parcel2source.pos = parcel_time.brainordinate.pos * 10;     %units: cm -> mm
parcel2source.inside = parcel_time.brainordinate.mask ~= 0; %inside brain positions
parcel2source.unit = 'mm';

%t-values: code group comparison
disp( strcat(num2str((1:length(groupnames))'),': ',groupnames) )
%use indices to refer to groups
CodeA = {[1 2],[1 3],[1 4],[5 4],[5 2],[5 3]};
%color
Cstat = cbrewer('div','RdBu',64);

%t-values computation & plot
for m = 1:length(CodeA)
    
    %group code
    code = CodeA{m};
    %t-value
    [~,~,~,tstat] = ttest2( CoM(:,indices.(groupnames{code(1)}))', CoM(:,indices.(groupnames{code(2)}))' );
    
    %t-value -> structure
    tval = tstat.tstat;
    
    %put CoM from parcel to sources
    parcel2source.tval = parcel_time.brainordinate.mask;
    for n = 1:length(tval)
        %exchange the mask index with the related frequency band power of that mask
        parcel2source.tval(parcel2source.tval == n) = tval(n);
    end
    
    %interpolate
    cfg = [];
    cfg.parameter = 'tval';
    cfg.downsample = 2;
    cfg.method = 'cubic';
    parcel2source_interpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
    
    %plot t-values
    cfg = [];
    cfg.method = 'surface';
    cfg.funparameter = 'tval';
    cfg.projmethod = 'nearest';
    cfg.funcolormap = Cstat;
    cfg.camlight = 'no';
    cfg.funcolorlim = [-4,4];
    ft_sourceplot(cfg,parcel2source_interpol)
    gr1_name = replace(groupnames{code(1)},'_','. ');
    gr2_name = replace(groupnames{code(2)},'_','. ');
    set(gcf,'Position',[10 10 650 700])
    title(['CoE: ',upper(gr1_name(1)),gr1_name(2:end),' vs. ',upper(gr2_name(1)),gr2_name(2:end)],'FontSize',32)
    colorbar('off')
    view([-90 0])
    c = colorbar('SouthOutside');
    c.FontSize = 32;
    c.FontWeight = 'bold';
    c.Label.String = {[upper(gr1_name),' < ',upper(gr2_name),'                   ',upper(gr2_name),' < ',upper(gr1_name)],...
        [''],...
        [''],...
        ['T-score']};
    c.Label.Position = [0 4 0];
    %left view
    print([mpath,'/parcel/ana_CoM/',replace(groupnames{code(1)},'_',' '),'VS',replace(groupnames{code(2)},'_',' '),'hz_left.tiff'],'-dtiff','-r300');
    %right view
    view([90 0])
    print([mpath,'/parcel/ana_CoM/',replace(groupnames{code(1)},'_',' '),'VS',replace(groupnames{code(2)},'_',' '),'hz_right.tiff'],'-dtiff','-r300');
    %change views and make images (right -> left -> top)
    view([0 90])
    print([mpath,'/parcel/ana_CoM/',replace(groupnames{code(1)},'_',' '),'VS',replace(groupnames{code(2)},'_',' '),'hz_top.tiff'],'-dtiff','-r300');
    
    close
            
end
