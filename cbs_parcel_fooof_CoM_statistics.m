%% Statistical Testing on Parcels of CoM [Cluster Permutation Test]

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

%load cbs project infos
load([mpath,'/cbs_info.mat']);  %cbs patients info

%load data 'CoM'
load([mpath,'/parcel/ana_CoM/CoM.mat']);CoE = tmp; clear tmp %CoM Parcel
%load fooof power for parcels
load([mpath,'/parcel/ana_power/fooof_spectra.mat']);

%load parcellation
load([mpath,'/parcellation.mat']);
%template grid (for later plotting)
template_grid.pos = parcel.pos; %units: m -> cm
template_grid.labels = parcel.masklabel;

%plot statistics
load([ft_path,'/template/anatomy/surface_white_both.mat'])
mesh = ft_convert_units(mesh,'cm');

%store figures
if ~exist([mpath,'/parcel/ana_stat/'])
    mkdir([mpath,'/parcel/ana_stat/'])
end

%define parcel neighbours
neigh = cbs_prepare_neighbours(parcel,1,'parcel','same');

%groupnames
groupnames = fieldnames(CoE.indices);

%color for statistical value
Cstat = cbrewer('div','RdBu',64);

%code group comparison (statistics)
disp( strcat(num2str((1:length(groupnames))'),': ',groupnames) )
%use indices to groups and compare means
CodeA = {[1 2],[1 3],[1 4],[1 5],[5 2],[5 3],[5 4],[2 3]};
%store eventually significant regions (each group comparison in one column)
SigR = table();

p_v = [];

for j = 1:length(CodeA)
    
    %group code
    code = CodeA{j};
    %group index
    gr1 = CoE.indices.( groupnames{code(1)} );
    gr2 = CoE.indices.( groupnames{code(2)} );
    %group CoE
    gr1_CoE = CoE.CoM(:,gr1);
    gr2_CoE = CoE.CoM(:,gr2);
    
    %labels
    labels = CoE.labels;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Statistical Test %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %hack into a fieldtrip-like data
    gr1_data = {};
    tmp = [];
    for k = 1:size(gr1_CoE,2)
        tmp.CoE = squeeze(gr1_CoE(:,k));
        tmp.label = labels;
        %act as if it was 'power'
        tmp.freq = 5;
        tmp.dimord = 'chan_freq';
        %put tmp in data
        gr1_data{k} = tmp;
    end
    %hack into a fieldtrip-like data
    gr2_data = {};
    tmp = [];
    for k = 1:size(gr2_CoE,2)
        tmp.CoE = squeeze(gr2_CoE(:,k));
        tmp.label = labels;
        %act as if it was 'power'
        tmp.freq = 5;
        tmp.dimord = 'chan_freq';
        %put tmp in data
        gr2_data{k} = tmp;
    end
    
    %___ Statistical Test [CoE]
    
    %group frequencies (for statistical testing)
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'CoE';
    gr1_data = ft_freqgrandaverage(cfg,gr1_data{:});
    gr2_data = ft_freqgrandaverage(cfg,gr2_data{:});
    
    %cluster permutation test
    cfg = [];
    cfg.neighbours = neigh;
    cfg.minnbchan = 1;
    cfg.method = 'montecarlo';
    cfg.statistic = 'ft_statfun_indepsamplesT';
    cfg.parameter = 'CoE';
    cfg.correctm = 'cluster';
    cfg.numrandomization = 10000;
    cfg.alpha = 0.025;
    cfg.clusteralpha = 0.025;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail = 0;
    %statistical parameter settings
    design(1,:) = [ones(1,size(gr1_data.CoE,1))*1,ones(1,size(gr2_data.CoE,1))*2];
    cfg.design = design;
    cfg.ivar = 1;
    %perform statistical test
    s = ft_freqstatistics(cfg,gr1_data,gr2_data);
    %clear 'design'
    clear design
    
    %store t-scores that are part of a 'significant' cluster
    SigR.([groupnames{code(1)},'VS',groupnames{code(2)}]) = s.stat .* s.mask;
    
    if sum(s.mask) ~= 0
        
        p_v(j) = min(s.prob);
        
        %extract only relevant t-values
        tval = s.stat .* s.mask;
        
        %make a tmp structure
        tmp.pos = template_grid.pos;
        
        %put tval from parcel to sources
        tmp.tval = parcel.mask;
        for n = 1:length(tval)
            %exchange the mask index with the related frequency band power of that mask
            tmp.tval(tmp.tval == n) = tval(n);
        end
        
        %interpolate
        cfg = [];
        cfg.parameter = 'tval';
        cfg.downsample = 2;
        cfg.method = 'cubic';
        tmp_interpol = ft_sourceinterpolate(cfg,tmp,mesh);
        
        %plot t-values
        cfg = [];
        cfg.method = 'surface';
        cfg.funparameter = 'tval';
        cfg.projmethod = 'nearest';
        cfg.funcolormap = Cstat;
        cfg.camlight = 'no';
        cfg.funcolorlim = [-3,3];
        cfg.colorbar = 'no';
        ft_sourceplot(cfg,tmp_interpol)
        set(gcf,'Position',[10 10 650 700])
        title(['CoE: ',upper(groupnames{code(1)}),' vs. ',upper(groupnames{code(2)})],'FontSize',32)
        colorbar('off')
        view([-90 0])
        c = colorbar('SouthOutside');
        c.FontSize = 32;
        c.FontWeight = 'bold';
        c.Label.String = {[upper(groupnames{code(1)}),' < ',upper(groupnames{code(2)}),'                   ',upper(groupnames{code(2)}),' < ',upper(groupnames{code(1)})],...
            [''],...
            [''],...
            ['T-score']};
        c.Label.Position = [0 4 0];
        %left view
        print([mpath,'/parcel/ana_stat/',groupnames{code(1)},'_',groupnames{code(2)},'_CoM_stat_left.tiff'],'-dtiff','-r300');
        %right view
        view([90 0])
        print([mpath,'/parcel/ana_stat/',groupnames{code(1)},'_',groupnames{code(2)},'_CoM_stat_right.tiff'],'-dtiff','-r300');
        %change views and make images (right -> left -> top)
        view([0 90])
        print([mpath,'/parcel/ana_stat/',groupnames{code(1)},'_',groupnames{code(2)},'_CoM_stat_top.tiff'],'-dtiff','-r300');

    else
        p_v(j) = nan;
        
        warning(['No significant clusters | ',replace(groupnames{code(1)},'_',' '),' vs ',replace(groupnames{code(2)},'_',' '),' | CoM'])
    end
    
end

save([mpath,'/parcel/ana_stat/Significant_Regions.mat'],'SigR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual CoE per parcel in ROI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use overlapping clusters of hc vs. atyp & pd vs. atyp
parcels_of_interest = parcel.masklabel( logical( SigR.hcVSaps .* SigR.pdVSaps ) );
%related CoE values
CoE_ROI = CoE.CoM( logical( SigR.hcVSaps .* SigR.pdVSaps ),: );
%seperate CoE into groups
CoE_ROI_hc = mean( CoE.CoM( logical( SigR.hcVSaps .* SigR.pdVSaps ),CoE.indices.hc ) ,1 )';
CoE_ROI_pd = mean( CoE.CoM( logical( SigR.hcVSaps .* SigR.pdVSaps ),CoE.indices.pd ) ,1 )';
CoE_ROI_atyp = mean( CoE.CoM( logical( SigR.hcVSaps .* SigR.pdVSaps ),CoE.indices.aps ) ,1 )';

jitter = -0.15:0.005:0.15;

hc_txt = CoE.sub(CoE.indices.hc);
pd_txt = CoE.sub(CoE.indices.pd);
atyp_txt = CoE.sub(CoE.indices.aps);

col = cbrewer('qual', 'Set3', 12, 'pchip');

figure('Position',[25 25 600 600])
h = boxplot(vertcat(CoE_ROI_hc,CoE_ROI_pd,CoE_ROI_atyp),...
            vertcat(repmat(1,sum(CoE.indices.hc),1),2*repmat(1,sum(CoE.indices.pd),1),3*repmat(1,sum(CoE.indices.aps),1)),...
            'Color',[col(1,:);col(6,:);col(4,:)],...
            'symbol','');
set(h,'linew',1)
hold on;
scatter(ones(1,length(CoE_ROI_hc)) + jitter(1,randperm(length(jitter),sum(CoE.indices.hc))),CoE_ROI_hc,50,'MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:)); hold on;
scatter(2*ones(1,length(CoE_ROI_pd)) + jitter(1,randperm(length(jitter),sum(CoE.indices.pd))),CoE_ROI_pd,50,'MarkerFaceColor',col(6,:),'MarkerEdgeColor',col(6,:)); hold on;
scatter(3*ones(1,length(CoE_ROI_atyp)) + jitter(1,randperm(length(jitter),sum(CoE.indices.aps))),CoE_ROI_atyp,50,'MarkerFaceColor',col(4,:),'MarkerEdgeColor',col(4,:)); hold on;
ylabel('CoE')
ylim([11.5,25])
drawbrace([1 24],[3 24],5,'Color',[0 0 0],'LineWidth',2); text([1.933 1.933],[24.5 24.5],'**','FontSize',20);
drawbrace([2 23],[3 23],5,'Color',[0 0 0],'LineWidth',2); text([2.433 2.433],[23.5 23.5],'**','FontSize',20);
set(gca,'XTickLabel','')
set(gca,'FontSize',20)

print(gcf,[mpath,'/parcel/ana_CoM/CoE_ROI.tiff'],'-dtiff','-r300')

close all

%compute statistic

%Age of groups
age_atyp = []; age_hc = []; age_pd = []; age_psp = []; age_cbs = [];
%catch subjects with problems
catch_subs = [];
subjects = CoE.sub( logical( CoE.indices.aps + CoE.indices.hc + CoE.indices.pd ) );
for i = 1:length(subjects)
    
    %for some Hc the Age was directly entered in the structure under '.age'
    if contains( 'age',fieldnames(cbs_info.(subjects{i})) )
        
        current_age = cbs_info.(subjects{i}).age;
        
        if isstr(current_age)
            current_age = str2num(current_age);
        end
        
    else
        
        %birthday and recording date
        birthday = cbs_info.(subjects{i}).birth;
        recording_day = cbs_info.(subjects{i}).date;
        
        if isempty(birthday)
            catch_subs = vertcat(catch_subs,subjects(i));
            continue
        end
        
        %get age
        current_age = cbs_GetAge(birthday,recording_day);
        
    end
    
    if ismember( subjects{i}, CoE.sub(CoE.indices.hc) )
       %put current_age
       age_hc = vertcat(age_hc,current_age);
    end
    
    if ismember( subjects{i}, CoE.sub(CoE.indices.pd) )
       %put current_age
       age_pd = vertcat(age_pd,current_age);
    end
    
    if ismember( subjects{i}, CoE.sub(CoE.indices.aps) )
       %put current_age
       age_atyp = vertcat(age_atyp,current_age);
    end
    
    if ismember( subjects{i}, CoE.sub(CoE.indices.cbs) )
       %put current_age
       age_cbs = vertcat(age_cbs,current_age);
    end

    if ismember( subjects{i}, CoE.sub(CoE.indices.psp) )
       %put current_age
       age_psp = vertcat(age_psp,current_age);
    end
    
end

age_ = vertcat(age_atyp,age_hc,age_pd);
age_ = age_ - mean(age_);

groups = vertcat(repmat({'atyp'},length(CoE_ROI_atyp),1),repmat({'hc'},length(CoE_ROI_hc),1),repmat({'pd'},length(CoE_ROI_pd),1));
groups = categorical(groups);

%make a table for linear model
C = table( vertcat(CoE_ROI_atyp,CoE_ROI_hc,CoE_ROI_pd),groups,age_,'VariableNames',{'CoE','group','age'} );

%linear model
mdl = fitlm(C,'CoE ~ 1 + group + age + group*age');

%save p-value
disp( mdl.Coefficients.pValue(2) );
disp( mdl.Coefficients.pValue(3) );

%mean CoE
disp(['Atyp CoE ROI | Mean: ',num2str(mean(CoE_ROI_atyp)),' | Std: ',num2str(std(CoE_ROI_atyp))])
disp(['HC CoE ROI | Mean: ',num2str(mean(CoE_ROI_hc)),' | Std: ',num2str(std(CoE_ROI_hc))])
disp(['PD CoE ROI | Mean: ',num2str(mean(CoE_ROI_pd)),' | Std: ',num2str(std(CoE_ROI_pd))])

%%%%%%%%%%%%%%%%%%%%%%
%%% ROI Topography %%%
%%%%%%%%%%%%%%%%%%%%%%

Cblues = cbrewer('seq','Blues',64,'cubic');

%interesting labels (contribute to clusters in statistic)
roi = find( abs(SigR.hcVSaps .* SigR.pdVSaps) > 0 );

%make a tmp structure
tmp.pos = template_grid.pos;

%put roi from parcel to sources
tmp.roi = 0.2*ones(length(parcel.mask),1);
for n = 1:length(roi)
    %exchange the mask index with the related frequency band power of that mask
    tmp.roi(parcel.mask == roi(n)) = 1;
end

%interpolate
cfg = [];
cfg.parameter = 'roi';
cfg.downsample = 2;
cfg.method = 'cubic';
tmp_interpol = ft_sourceinterpolate(cfg,tmp,mesh);

%plot t-values
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'roi';
cfg.projmethod = 'nearest';
cfg.funcolormap = Cblues;
cfg.camlight = 'no';
cfg.funcolorlim = [0,1];
cfg.colorbar = 'no';
ft_sourceplot(cfg,tmp_interpol)
set(gcf,'Position',[10 10 650 700])
%print it
print([mpath,'/parcel/ana_stat/ROI_cluster_top.tiff'],'-dtiff','-r300');
