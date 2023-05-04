%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Group Spectral Feature analysis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clean workspace and command window
clearvars;
clc;

%%%%%%%%%%%%%%%%%%%%
%%%___settings___%%%
%%%%%%%%%%%%%%%%%%%%

%path settings
mpath = 'C:/data';                                              %mainpath
ft_path = 'C:/toolboxes/fieldtrip-20201214';                    %fieltrip path
fct_path = [mpath,'/functions'];                                %function path (my own functions)
scp_path = [mpath,'/scripts'];                                  %script path
plt_path = 'C:/toolboxes/RainCloudPlots-master/tutorial_matlab';%plot path

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path,plt_path);
ft_defaults;

%load cbs project infos, cortex-areal sensor infos, define condition
load([mpath,'/cbs_info.mat']);  %cbs patients info

%load source information (for positions, unit, labels)
load([mpath,'/parcellation.mat']);

%load Spectral Features (Peak Oscillation & Amplitude)
load([mpath,'/parcel/ana_shift/spectral_features.mat']);

%load significant CoE clusters
load([mpath,'/parcel/ana_stat/Significant_Regions.mat']);

%subjects
subjects = SpecFeat.subjects;

%tremor subjects (from pd)
tremor = cbs_clean_subjects(subjects,cbs_info,{'rest','tremor'},'yes');
%notremor subjects
notremor = cbs_clean_subjects(subjects,cbs_info,{'rest','tremor'},'-yes');

%%%%%%%%%%%%%%
%%% Script %%%
%%%%%%%%%%%%%%

%groups
groups = {'hc',{'cbs','psp'},'pd'};

%group indices
indices.hc = contains(subjects,'hc');
indices.aps = contains(subjects,{'cbs','psp'});
indices.ips = logical( contains(subjects,'pd') .* contains(subjects,notremor) );

%use overlapping clusters of hc vs. atyp & pd vs. atyp
parcels_of_interest = parcel.masklabel( logical( SigR.hcVSaps .* SigR.pdVSaps ) );

%choose major oscillation by amplitude per parcel and per subject
use_amp = 'no'; %use 'no' & 'yes'
%make videos
use_video = 'yes';
%color
col = cbrewer('qual', 'Set3', 12, 'pchip');

%choose frequency bin interval for oscillations
intervall = [4,30];

%some variables that are used for plotting & analysing
Osc_hc = []; Osc_atyp = []; Osc_pd = [];                                % <- Get Center Frequencies of Oscillations of the group (cumulative over all parcels of interest)
Amp_hc = []; Amp_atyp = []; Amp_pd = [];                                % <- Get absolute Power / Amplitude of the group (cumulative over all parcels of interest)

%get all differences as t-value pwer parcel (but ... mostly to fews osillations in a single parcel to extract some meaning here)
Tval_Osc = zeros(3,length(parcels_of_interest));
Tval_Amp = zeros(3,length(parcels_of_interest));

%Oscillations per group and parcel
GroupOsc.Osc = [];
GroupOsc.rows = {'hc','aps','ips'}';
GroupOsc.parcels = parcels_of_interest;
GroupOsc.indices = indices;

for a = 1:length(parcels_of_interest)
    %area
    area  = find(contains(SpecFeat.labels,parcels_of_interest{a}));
    %oscillation data
    Os = cell(length(groups),1);
    %amplitude data
    Am = cell(length(groups),1);
    
    for n = 1:length(groups)
        
        %choose group
        gr = groups{n};
        
        for k = 1:length(subjects)
            if contains(subjects(k),gr) && ~contains(cbs_info.(subjects{k}).rest.tremor,'yes')
                
                %oscillations in 'area' for subject k
                osc_candidates = SpecFeat.Osc{ ismember(SpecFeat.subjects,subjects{k}) }{area};
                
                %amplitude to the oscillations in 'area' for subject k
                amp_candidates = SpecFeat.Amp{ ismember(SpecFeat.subjects,subjects{k}) }{area};
                
                if length(intervall) == 2
                    %apply intervall restrictions
                    idx = logical( ( osc_candidates > intervall(1) ) .* ( osc_candidates < intervall(2) ) );
                    %prune osc_candidates
                    osc_candidates = osc_candidates(idx);
                    %prune amp_candidates
                    amp_candidates = amp_candidates(idx);
                end
                
                if strcmpi(use_amp,'yes')
                    %select the major ('maximal') oscillation
                    [amp_candidates,ch] = max(amp_candidates);
                    %choose oscillation
                    osc_candidates = osc_candidates(ch);
                end
                
                if isrow(osc_candidates)
                    osc_candidates = osc_candidates';
                end
                if isrow(amp_candidates)
                    amp_candidates = amp_candidates';
                end                    
                    
                %save oscillations and amplitudes
                Os{n} = vertcat(Os{n},osc_candidates);
                Am{n} = vertcat(Am{n},amp_candidates);
                
            end
        end
        
    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Oscillatory Center Frequency %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Center Frequency of Oscillations
        Osc_hc  = vertcat( Osc_hc,Os{1} );
        Osc_atyp= vertcat( Osc_atyp,Os{2} );
        Osc_pd  = vertcat( Osc_pd,Os{3} );
        
        %store for parcel oscillations for statistical testing
        GroupOsc.Osc = horzcat(GroupOsc.Osc,Os);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Amplitudes in defined Intervall %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Amplitude of Center Frequency
        Amp_hc  = vertcat( Amp_hc,Am{1} );
        Amp_atyp= vertcat( Amp_atyp,Am{2} );
        Amp_pd  = vertcat( Amp_pd,Am{3} );
        
end

%save structure
if all(ismember(parcel.masklabel,parcels_of_interest))
    save([mpath,'/parcel/ana_shift/GroupOsc.mat'],'GroupOsc');
else
    save([mpath,'/parcel/ana_shift/GroupOscCluster.mat'],'GroupOsc');
end
clear GroupOsc

%just to adjust the Raincloud Plot
if intervall(1) == 4 && intervall(2) == 30 && strcmpi(use_amp,'no')
    yli = [-0.05 0.075; -4.25 7];
elseif intervall(1) == 4 && intervall(2) == 30 && strcmpi(use_amp,'yes')
    yli = [-0.07 0.1; -4 5.5];
elseif intervall(1) == 4 && intervall(2) == 13 && strcmpi(use_amp,'no')
    yli = [-0.15 0.25; -6 7.5];
elseif intervall(1) == 4 && intervall(2) == 13 && strcmpi(use_amp,'yes')
    yli = [-0.15 0.25; -5 7.5];
elseif intervall(1) == 13 && intervall(2) == 30 && strcmpi(use_amp,'no')
    yli = [-0.07 0.12; -3 4.5];
elseif intervall(1) == 13 && intervall(2) == 30 && strcmpi(use_amp,'yes')
    yli = [-0.1 0.15; -2.5 3.75];
end

%Raincloud Plots
figure('Position',[50 50 1300 400])
subplot(1,2,1)
h1 = raincloud_plot(Osc_hc, 'box_on', 1, 'color', col(1,:), 'alpha', 0.2,...
    'box_dodge', 1, 'box_dodge_amount', 0.15, 'dot_dodge_amount', .15,'box_col_match', 0 );
h1{1}.EdgeColor = col(1,:);
h1{3}.EdgeColor = col(1,:);
h1{4}.Color = col(1,:);
h1{5}.Color = col(1,:);
h1{6}.Color = col(1,:);
h1{2}.MarkerFaceAlpha = 0.33;
h2 = raincloud_plot(Osc_atyp, 'box_on', 1, 'color', col(4,:), 'alpha', 0.2,...
    'box_dodge', 1, 'box_dodge_amount', 0.35, 'dot_dodge_amount', .35,'box_col_match', 0 );
h2{1}.EdgeColor = col(4,:);
h2{1}.EdgeColor = col(4,:);
h2{3}.EdgeColor = col(4,:);
h2{4}.Color = col(4,:);
h2{5}.Color = col(4,:);
h2{6}.Color = col(4,:);
h2{2}.MarkerFaceAlpha = 0.33;
h3 = raincloud_plot(Osc_pd, 'box_on', 1, 'color', col(6,:), 'alpha', 0.2,...
    'box_dodge', 1, 'box_dodge_amount', 0.55, 'dot_dodge_amount', .55,'box_col_match', 0 );
h3{1}.EdgeColor = col(6,:);
h3{1}.EdgeColor = col(6,:);
h3{3}.EdgeColor = col(6,:);
h3{4}.Color = col(6,:);
h3{5}.Color = col(6,:);
h3{6}.Color = col(6,:);
h3{2}.MarkerFaceAlpha = 0.33;
%legend([h1{1} h2{1} h3{1}], {'Hc', 'Atyp. PD', 'Idio. PD'})
title('Peak frequency','FontSize',22,'Fontweight','bold')
if intervall(1) == 4 && intervall(2) == 30
    set(gca,'XLim', [2 32], 'YLim', yli(1,:), 'FontSize',20);
elseif intervall(1) == 4 && intervall(2) == 13
    set(gca,'XLim', [3 14], 'YLim',yli(1,:), 'FontSize',20);
elseif intervall(1) == 13 && intervall(2) == 30
    set(gca,'XLim', [11 33], 'YLim',yli(1,:), 'FontSize',20);
end
set(gca,'XTick')
xlabel('Frequency','FontSize',20)
ylabel('Probabilidy Density','FontSize',20)
box off

%Inbetween catch the 'sink' of the atyp. group
LocalMin_atyp = round( ( h2{1}.XData( islocalmin(h2{1}.YData) ) ) ,1);
%Inbetween catch the 'sink' of the atyp. group
LocalMin_hc = round( ( h1{1}.XData( islocalmin(h1{1}.YData) ) ) ,1);
%Inbetween catch the 'sink' of the atyp. group
LocalMin_pd = round( ( h3{1}.XData( islocalmin(h3{1}.YData) ) ) ,1);

%Inbetween catch the 'sink' of the atyp. group
LocalMax_atyp = round( ( h2{1}.XData( islocalmax(h2{1}.YData) ) ) ,1);
%Inbetween catch the 'sink' of the atyp. group
LocalMax_hc = round( ( h1{1}.XData( islocalmax(h1{1}.YData) ) ) ,1);
%Inbetween catch the 'sink' of the atyp. group
LocalMax_pd = round( ( h3{1}.XData( islocalmax(h3{1}.YData) ) ) ,1);

subplot(1,2,2)
h1 = raincloud_plot(Amp_hc, 'box_on', 1, 'color', col(1,:), 'alpha', 0.2,...
    'box_dodge', 1, 'box_dodge_amount', 0.15, 'dot_dodge_amount', .15,'box_col_match', 0 );
h1{1}.EdgeColor = col(1,:);
h1{3}.EdgeColor = col(1,:);
h1{4}.Color = col(1,:);
h1{5}.Color = col(1,:);
h1{6}.Color = col(1,:);
h1{2}.MarkerFaceAlpha = 0.33;
h2 = raincloud_plot(Amp_atyp, 'box_on', 1, 'color', col(4,:), 'alpha', 0.2,...
    'box_dodge', 1, 'box_dodge_amount', 0.35, 'dot_dodge_amount', .35,'box_col_match', 0 );
h2{1}.EdgeColor = col(4,:);
h2{1}.EdgeColor = col(4,:);
h2{3}.EdgeColor = col(4,:);
h2{4}.Color = col(4,:);
h2{5}.Color = col(4,:);
h2{6}.Color = col(4,:);
h2{2}.MarkerFaceAlpha = 0.33;
h3 = raincloud_plot(Amp_pd, 'box_on', 1, 'color', col(6,:), 'alpha', 0.2,...
    'box_dodge', 1, 'box_dodge_amount', 0.55, 'dot_dodge_amount', .55,'box_col_match', 0 );
h3{1}.EdgeColor = col(6,:);
h3{1}.EdgeColor = col(6,:);
h3{3}.EdgeColor = col(6,:);
h3{4}.Color = col(6,:);
h3{5}.Color = col(6,:);
h3{6}.Color = col(6,:);
h3{2}.MarkerFaceAlpha = 0.33;
lgd = legend([h1{1} h2{1} h3{1}], {'HC', 'APS', 'PD'});

title('Peak amplitude','FontSize',22,'Fontweight','bold')
if intervall(1) == 4 && intervall(2) == 30
    set(gca,'XLim', [-0.05 0.9], 'YLim',yli(2,:), 'FontSize',20);
elseif intervall(1) == 4 && intervall(2) == 13
    set(gca,'XLim', [-0.1 0.9], 'YLim',yli(2,:), 'FontSize',20);
elseif intervall(1) == 13 && intervall(2) == 30
    set(gca,'XLim', [-0.1 0.9], 'YLim',yli(2,:), 'FontSize',20);
end
set(gca,'XTick')
xlabel('Amplitude','FontSize',20)
ylabel('Probabilidy Density','FontSize',20)
box off

if intervall(1) == 4 && intervall(2) == 30
   print([mpath,'/parcel/ana_shift/Oscillations_AmplitudeUse_',use_amp,'.tiff'],'-dtiff','-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ---------------------------- Typ I ---------------------------------%%%
%%% Regional Approach | y-Axis Gradient: Center Frequency and Amplitude %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tmp

%load mesh
load([ft_path,'/template/anatomy/surface_white_both.mat'])
mesh = ft_convert_units(mesh,'cm');
%load cortical grid
load([mpath,'\cortical_grid.mat'])
template_grid.pos = cortical_grid' .* 100; %units: m -> cm

%Parcel Coordinates
ParcelCoords = cell(length(parcel.masklabel),1);
%Parcel Centroid
ParcelCentroid = zeros(length(ParcelCoords),3);
for p = 1:size(ParcelCoords,1)
    ParcelCoords{p} = parcel.pos( parcel.mask == p , : );
    ParcelCentroid(p,:) = mean( ParcelCoords{p} );
end

%tmp
tmp.pos = parcel.pos;
tmp.unit = parcel.unit;

%color
col = cbrewer('qual', 'Set3', 12, 'pchip');

%%% --- Group variables as a function of y-coordinates --- %%%

%Abs_Oscillations & Standard Error
Abs_avg_All = []; Abs_All_IQR = [];
%Median(Oscillations) & .25 / .75-Qauntiles
OscMed = []; OscIQR= []; OscAvg = []; OscStd = [];
%Average(Amplitudes) & Standard Error
AmpMed = []; AmpAvg = []; AmpIQR = [];
%All Cluster Coordinates used
AllClusterCoords = {};
%All Parcels used
AllParcelsOfInterest = {};

%Oscillations per region and group
OscRegionHc = [];
OscRegionAtyp = [];
OscRegionPd = [];

xdata_hc = [];
xdata_atyp = [];
xdata_pd = [];
ydata_hc = [];
ydata_atyp = [];
ydata_pd = [];

%y-coordinates of the cortical_grid
ycoords = unique(template_grid.pos(:,2));
%delete first y-coordinate as it is only 2 mm away from the following one
ycoords = ycoords(2:end);

%margin to include region if parcel centroid is within
y_margin = 2.5;

%exclude parcels with certain label
exclude_parcels = {'Cerebellum'};
%already compute their corresponding mask number
idx_exclude = find(contains(parcel.masklabel,exclude_parcels));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EITHER RUN THIS OR ... %%% NOsc, Peak Frequency, Peak Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([mpath,'/parcel/time/hc01/hc01_parcel_time.mat']);
parcel2source.pos = parcel_time.brainordinate.pos * 10;     %units: cm -> mm
parcel2source.inside = parcel_time.brainordinate.mask ~= 0; %inside brain positions
parcel2source.mask = parcel_time.brainordinate.mask;        %mask
parcel2source.unit = 'mm';

if strcmpi(use_amp,'no')
    figure('Position',[50 50 1200 600])
    %group parcels together [force if _L ... appear then also _R and vice versa]
    for par = 1 : length(ycoords)
        
        %store parcels for cluster
        clu = [];
        %apply yPos to Parcel Centroids
        for c = 1:length(ParcelCentroid)
            
            yParcel = unique( ParcelCentroid(c,2) );
            
            %distance (cm) to any of the choose y positions
            distances = abs( yParcel - ycoords(par) ) <= y_margin;
            %check if all three coordinates are in the parcel
            check = sum(sum(distances,2),1) >= min( [3,length(yParcel)] );
            
            if any(check)
                %if a parcel fullfills the 'y - criteria' store its index
                clu = cat(1,clu,c);
                %and take its friend on the other hemisphere
                if endsWith( parcel.masklabel{c},'R')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_R','_L')) ) );
                    clu = cat(1,clu,friend_index);
                elseif endsWith( parcel.masklabel{c},'L')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_L','_R')) ) );
                    clu = cat(1,clu,friend_index);
                end
            end
            
        end
        
        %exclude parcels ...
        clu = unique( clu(~ismember(clu,idx_exclude)) );
        %and in the last-1 step also disregard Frontal_Sup. Looks MUCH better then
        if par == length(ycoords)-1; idx_frontal_sup = [find(contains(parcel.masklabel,'Frontal_Sup_L')),find(contains(parcel.masklabel,'Frontal_Sup_R'))]; clu = clu(~ismember(clu,idx_frontal_sup)) ;end
        
        %Parcel Labels
        parcel2source.reg = ismember(parcel2source.mask,clu);
        
        %highlight region 'reg' on mesh -> on mesh
        cfg = [];
        cfg.parameter = 'reg';
        cfg.downsample = 1;
        cfg.method = 'cubic';
        SourceInterpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
        %highlight the positions for the specific region
        highlight.pos = SourceInterpol.pos(logical(SourceInterpol.reg),:);
        
        h1 = subplot(3,2,[1 3]);
        %basic background plot
        ft_plot_mesh( mesh, 'facealpha',0.1 ); hold on;
        %add region in another color
        ft_plot_mesh( highlight, 'vertexcolor',[.8 .8 .8] )
        view([90,0])
        
        %current parcels of interest
        parcels_of_interest = parcel.masklabel(clu);
        
        %some variables that are used for plotting & analysing
        Osc_hc = []; Osc_atyp = []; Osc_pd = [];                                % <- Get Center Frequencies of Oscillations of the group (cumulative over all parcels of interest)
        Amp_hc = []; Amp_atyp = []; Amp_pd = [];                                % <- Get absolute Power / Amplitude of the group (cumulative over all parcels of interest)
        Abs_hc = []; Abs_atyp = []; Abs_pd = [];                                % <- Get absolute number of Oscillations of the group (Number of Oscillations per parcel and Group)
        
        for a = 1:length(parcels_of_interest)
            %area
            area  = find(contains(SpecFeat.labels,parcels_of_interest{a}));
            %oscillation data
            Os = cell(length(groups),1);
            %amplitude data
            Am = cell(length(groups),1);
            
            for n = 1:length(groups)
                
                %choose group
                gr = groups{n};
                
                for k = 1:length(subjects)
                    if contains(subjects(k),gr) && ~contains(cbs_info.(subjects{k}).rest.tremor,'yes')
                        
                        %oscillations in 'area' for subject k
                        osc_candidates = SpecFeat.Osc{k}{area};
                        
                        %amplitude to the oscillations in 'area' for subject k
                        amp_candidates = SpecFeat.Amp{k}{area};
                        
                        if length(intervall) == 2
                            %apply intervall restrictions
                            idx = logical( ( osc_candidates > intervall(1) ) .* ( osc_candidates < intervall(2) ) );
                            %prune osc_candidates
                            osc_candidates = osc_candidates(idx);
                            %prune amp_candidates
                            amp_candidates = amp_candidates(idx);
                        end
                        
                        if strcmpi(use_amp,'yes')
                            %select the major ('maximal') oscillation
                            [amp_candidates,ch] = max(amp_candidates);
                            %choose oscillation
                            osc_candidates = osc_candidates(ch);
                        end
                        
                        if isrow(osc_candidates)
                            osc_candidates = osc_candidates';
                        end
                        if isrow(amp_candidates)
                            amp_candidates = amp_candidates';
                        end
                        
                        %save oscillations and amplitudes
                        Os{n} = vertcat(Os{n},osc_candidates);
                        Am{n} = vertcat(Am{n},amp_candidates);
                        
                    end
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Oscillatory Center Frequency %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Absolute number of oscillations in this area and intervall for the different groups
            Abs_hc(a) = length(Os{1});
            Abs_atyp(a) = length(Os{2});
            Abs_pd(a) = length(Os{3});
            
            %Center Frequency of Oscillations
            Osc_hc  = vertcat( Osc_hc,Os{1} );
            Osc_atyp= vertcat( Osc_atyp,Os{2} );
            Osc_pd  = vertcat( Osc_pd,Os{3} );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Amplitudes in defined Intervall %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %get Power for Parcel 'a' for all three groups as three elements of the array Pow
            Pow = cellfun(@sum,Am);
            
            %Amplitude of Center Frequency
            Amp_hc  = vertcat( Amp_hc,Am{1} );
            Amp_atyp= vertcat( Amp_atyp,Am{2} );
            Amp_pd  = vertcat( Amp_pd,Am{3} );
            
            %Average Amplitude per parcel and subject as given by group identification
            Avg_Pow_hc  = sum(Amp_hc) / (a * sum(indices.hc));
            Avg_Pow_atyp= sum(Amp_atyp) / (a * sum(indices.aps));
            Avg_Pow_pd  = sum(Amp_pd) /(a * sum(indices.ips));
            
            %clear Pow for security reasons ... ;)
            clear Pow
            
        end
        
        %Average Number of Oscillations per parcel and subject as given by group identification
        Avg_N_Osc_hc  = sum(Abs_hc)/ (length(parcels_of_interest) * sum(indices.hc));
        Avg_N_Osc_atyp = sum(Abs_atyp)/ (length(parcels_of_interest) * sum(indices.aps));
        Avg_N_Osc_pd  = sum(Abs_pd)/(length(parcels_of_interest) * sum(indices.ips));
        
        %store all Oscillations of Region a
        OscRegionHc{par}= Osc_hc;
        OscRegionAtyp{par}= Osc_atyp;
        OscRegionPd{par}= Osc_pd;
        
        %store Group Median Oscillations (for later Plotting)
        OscMed = vertcat(OscMed,[median(Osc_hc),median(Osc_atyp),median(Osc_pd)]);
        OscIQR = vertcat(OscIQR,[quantile(Osc_hc,.25),quantile(Osc_hc,.75),quantile(Osc_atyp,.25),quantile(Osc_atyp,.75),quantile(Osc_pd,.25),quantile(Osc_pd,.75)]);
        %store Group Mean Oscillations (for later Plotting)
        OscAvg = vertcat(OscAvg,[mean(Osc_hc),mean(Osc_atyp),mean(Osc_pd)]);
        OscStd = vertcat(OscStd,[std(Osc_hc),std(Osc_atyp),std(Osc_pd)]);
        
        %store Group Amplitude (for later Plotting)
        AmpMed = vertcat(AmpMed,[median(Amp_hc),median(Amp_atyp),median(Amp_pd)]);
        AmpAvg = vertcat(AmpAvg,[mean(Amp_hc),mean(Amp_atyp),mean(Amp_pd)]);
        AmpIQR = vertcat(AmpIQR,[quantile(Amp_hc,.25),quantile(Amp_hc,.75),quantile(Amp_atyp,.25),quantile(Amp_atyp,.75),quantile(Amp_pd,.25),quantile(Amp_pd,.75)]);
        %Absolute number of oscillations over y-axis
        Abs_avg_All = vertcat(Abs_avg_All,[Avg_N_Osc_hc,Avg_N_Osc_atyp,Avg_N_Osc_pd]);
        
        subplot(3,2,5)
        plot(ycoords(1:par),Abs_avg_All(:,1),'Color',col(1,:),'LineWidth',3)
        hold on
        plot(ycoords(1:par),Abs_avg_All(:,2),'Color',col(4,:),'LineWidth',3)
        hold on
        plot(ycoords(1:par),Abs_avg_All(:,3),'Color',col(6,:),'LineWidth',3)
        xlim([min(ycoords),max(ycoords)])
        ylim([2,2.5])
        xticklabels('')
        xlabel('y-coordinates')
        ylabel('Avg. N Oscillations')
        
        title('N Oscillations')
        
        h2 = subplot(3,2,[2 4]);
        r1 = raincloud_plot(Osc_hc, 'box_on', 1, 'color', col(1,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.12, 'dot_dodge_amount', .12,'box_col_match', 0 );
        r2 = raincloud_plot(Osc_atyp, 'box_on', 1, 'color', col(4,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.32, 'dot_dodge_amount', .32,'box_col_match', 0 );
        r3 = raincloud_plot(Osc_pd, 'box_on', 1, 'color', col(6,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.52, 'dot_dodge_amount', .52,'box_col_match', 0 );
        legend([r1{1} r2{1} r3{1}], {'Hc', 'Atyp. PD', 'Idio. PD'})
        title('Center of Frequency')
        set(gca,'XLim', [1,35]);
        %set(gca,'YLim', [-0.07,0.1]);
        xlabel('Frequency')
        ylabel('Probabilidy Density')
        box off
        
        %save density estimation
        xdata_hc{par} = r1{1}.XData; ydata_hc{par} = r1{1}.YData;
        xdata_atyp{par} = r2{1}.XData; ydata_atyp{par} = r2{1}.YData;
        xdata_pd{par} = r3{1}.XData; ydata_pd{par} = r3{1}.YData;
        
        %some jitter for groups & Amplitudes
        jit = -.21 : 0.01 : +.2;
        
        jit_hc_amp   = jit(randi([1,numel(jit)],length(Amp_hc),1))';
        jit_atyp_amp = jit(randi([1,numel(jit)],length(Amp_atyp),1))';
        jit_pd_amp   = jit(randi([1,numel(jit)],length(Amp_pd),1))';
        
        h3 = subplot(3,2,6);
        scatter(Amp_hc, 1 + jit_hc_amp,'MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'MarkerFaceAlpha',.5);
        hold on
        scatter(Amp_atyp, 2 + jit_atyp_amp,'MarkerEdgeColor',col(4,:),'MarkerFaceColor',col(4,:),'MarkerFaceAlpha',.5);
        hold on
        scatter(Amp_pd, 3 + jit_pd_amp,'MarkerEdgeColor',col(6,:),'MarkerFaceColor',col(6,:),'MarkerFaceAlpha',.5);
        line([mean(Amp_hc),mean(Amp_hc)],[ 1-0.2, 1+0.2],'Color',[0 0 0],'LineWidth',3)
        line([mean(Amp_atyp),mean(Amp_atyp)],[ 2-0.2, 2+0.2],'Color',[0 0 0],'LineWidth',3)
        line([mean(Amp_pd),mean(Amp_pd)],[ 3-0.2, 3+0.2],'Color',[0 0 0],'LineWidth',3)
        ylim([0,4])
        hold on
        
        ylabel('Groups')
        xlabel('Amplitude')
        
        yticklabels(' ')
        
        title('Amplitude of Oscillations')
        
        legend('Hc','Atyp. PD','Idio. PD')
        
        pause(1)
        
        if strcmpi(use_video,'yes')
            %get frame
            F{par} = getframe(gcf);
        end
        
        if par  < length(ycoords)
            cla(h1)
            cla(h2)
            cla(h3)
        end
        
    end
elseif strcmpi(use_amp,'yes')
    figure('Position',[50 50 1200 600])
    %group parcels together [force if _L ... appear then also _R and vice versa]
    for par = 1 : length(ycoords)
        
        %store parcels for cluster
        clu = [];
        %apply yPos to Parcel Centroids
        for c = 1:length(ParcelCentroid)
            
            yParcel = unique( ParcelCentroid(c,2) );
            
            %distance (cm) to any of the choose y positions
            distances = abs( yParcel - ycoords(par) ) <= y_margin;
            %check if all three coordinates are in the parcel
            check = sum(sum(distances,2),1) >= min( [3,length(yParcel)] );
            
            if any(check)
                %if a parcel fullfills the 'y - criteria' store its index
                clu = cat(1,clu,c);
                %and take its friend on the other hemisphere
                if endsWith( parcel.masklabel{c},'R')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_R','_L')) ) );
                    clu = cat(1,clu,friend_index);
                elseif endsWith( parcel.masklabel{c},'L')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_L','_R')) ) );
                    clu = cat(1,clu,friend_index);
                end
            end
            
        end
        
        %exclude parcels ...
        clu = unique( clu(~ismember(clu,idx_exclude)) );
        %and in the last-1 step also disregard Frontal_Sup. Looks MUCH better then
        if par == length(ycoords)-1; idx_frontal_sup = [find(contains(parcel.masklabel,'Frontal_Sup_L')),find(contains(parcel.masklabel,'Frontal_Sup_R'))]; clu = clu(~ismember(clu,idx_frontal_sup)) ;end
        
        %Parcel Labels
        parcel2source.reg = ismember(parcel2source.mask,clu);
        
        %highlight region 'reg' on mesh -> on mesh
        cfg = [];
        cfg.parameter = 'reg';
        cfg.downsample = 1;
        cfg.method = 'cubic';
        SourceInterpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
        %highlight the positions for the specific region
        highlight.pos = SourceInterpol.pos(logical(SourceInterpol.reg),:);
        
        h1 = subplot(3,2,[1 3 5]);
        %basic background plot
        ft_plot_mesh( mesh, 'facealpha',0.1 ); hold on;
        %add region in another color
        ft_plot_mesh( highlight, 'vertexcolor',[.8 .8 .8] )
        view([90,0])
        
        %current parcels of interest
        parcels_of_interest = parcel.masklabel(clu);
        
        %some variables that are used for plotting & analysing
        Osc_hc = []; Osc_atyp = []; Osc_pd = [];                                % <- Get Center Frequencies of Oscillations of the group (cumulative over all parcels of interest)
        Amp_hc = []; Amp_atyp = []; Amp_pd = [];                                % <- Get absolute Power / Amplitude of the group (cumulative over all parcels of interest)
        Abs_hc = []; Abs_atyp = []; Abs_pd = [];                                % <- Get absolute number of Oscillations of the group (Number of Oscillations per parcel and Group)
        
        for a = 1:length(parcels_of_interest)
            %area
            area  = find(contains(SpecFeat.labels,parcels_of_interest{a}));
            %oscillation data
            Os = cell(length(groups),1);
            %amplitude data
            Am = cell(length(groups),1);
            
            for n = 1:length(groups)
                
                %choose group
                gr = groups{n};
                
                for k = 1:length(subjects)
                    if contains(subjects(k),gr) && ~contains(cbs_info.(subjects{k}).rest.tremor,'yes')
                        
                        %oscillations in 'area' for subject k
                        osc_candidates = SpecFeat.Osc{k}{area};
                        
                        %amplitude to the oscillations in 'area' for subject k
                        amp_candidates = SpecFeat.Amp{k}{area};
                        
                        if length(intervall) == 2
                            %apply intervall restrictions
                            idx = logical( ( osc_candidates > intervall(1) ) .* ( osc_candidates < intervall(2) ) );
                            %prune osc_candidates
                            osc_candidates = osc_candidates(idx);
                            %prune amp_candidates
                            amp_candidates = amp_candidates(idx);
                        end
                        
                        if strcmpi(use_amp,'yes')
                            %select the major ('maximal') oscillation
                            [amp_candidates,ch] = max(amp_candidates);
                            %choose oscillation
                            osc_candidates = osc_candidates(ch);
                        end
                        
                        if isrow(osc_candidates)
                            osc_candidates = osc_candidates';
                        end
                        if isrow(amp_candidates)
                            amp_candidates = amp_candidates';
                        end
                        
                        %save oscillations and amplitudes
                        Os{n} = vertcat(Os{n},osc_candidates);
                        Am{n} = vertcat(Am{n},amp_candidates);
                        
                    end
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Oscillatory Center Frequency %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Absolute number of oscillations in this area and intervall for the different groups
            Abs_hc(a) = length(Os{1});
            Abs_atyp(a) = length(Os{2});
            Abs_pd(a) = length(Os{3});
            
            %Center Frequency of Oscillations
            Osc_hc  = vertcat( Osc_hc,Os{1} );
            Osc_atyp= vertcat( Osc_atyp,Os{2} );
            Osc_pd  = vertcat( Osc_pd,Os{3} );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Amplitudes in defined Intervall %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %get Power for Parcel 'a' for all three groups as three elements of the array Pow
            Pow = cellfun(@sum,Am);
            
            %Amplitude of Center Frequency
            Amp_hc  = vertcat( Amp_hc,Am{1} );
            Amp_atyp= vertcat( Amp_atyp,Am{2} );
            Amp_pd  = vertcat( Amp_pd,Am{3} );
            
            %Average Amplitude per parcel and subject as given by group identification
            Avg_Pow_hc  = sum(Amp_hc) / (a * sum(indices.hc));
            Avg_Pow_atyp= sum(Amp_atyp) / (a * sum(indices.aps));
            Avg_Pow_pd  = sum(Amp_pd) /(a * sum(indices.ips));
            
            %clear Pow for security reasons ... ;)
            clear Pow
            
        end
                
        %store all Oscillations of Region a
        OscRegionHc{par}= Osc_hc;
        OscRegionAtyp{par}= Osc_atyp;
        OscRegionPd{par}= Osc_pd;
        
        %store Group Median Oscillations (for later Plotting)
        OscMed = vertcat(OscMed,[median(Osc_hc),median(Osc_atyp),median(Osc_pd)]);
        OscIQR = vertcat(OscIQR,[quantile(Osc_hc,.25),quantile(Osc_hc,.75),quantile(Osc_atyp,.25),quantile(Osc_atyp,.75),quantile(Osc_pd,.25),quantile(Osc_pd,.75)]);
        %store Group Mean Oscillations (for later Plotting)
        OscAvg = vertcat(OscAvg,[mean(Osc_hc),mean(Osc_atyp),mean(Osc_pd)]);
        OscStd = vertcat(OscStd,[std(Osc_hc),std(Osc_atyp),std(Osc_pd)]);
        
        %store Group Amplitude (for later Plotting)
        AmpMed = vertcat(AmpMed,[median(Amp_hc),median(Amp_atyp),median(Amp_pd)]);
        AmpAvg = vertcat(AmpAvg,[mean(Amp_hc),mean(Amp_atyp),mean(Amp_pd)]);
        AmpIQR = vertcat(AmpIQR,[quantile(Amp_hc,.25),quantile(Amp_hc,.75),quantile(Amp_atyp,.25),quantile(Amp_atyp,.75),quantile(Amp_pd,.25),quantile(Amp_pd,.75)]);
        
        h2 = subplot(3,2,[2 4]);
        r1 = raincloud_plot(Osc_hc, 'box_on', 1, 'color', col(1,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.12, 'dot_dodge_amount', .12,'box_col_match', 0 );
        r2 = raincloud_plot(Osc_atyp, 'box_on', 1, 'color', col(4,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.32, 'dot_dodge_amount', .32,'box_col_match', 0 );
        r3 = raincloud_plot(Osc_pd, 'box_on', 1, 'color', col(6,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.52, 'dot_dodge_amount', .52,'box_col_match', 0 );
        legend([r1{1} r2{1} r3{1}], {'Hc', 'Atyp. PD', 'Idio. PD'})
        title('Center of Frequency')
        set(gca,'XLim', [1,35]);
        %set(gca,'YLim', [-0.07,0.1]);
        xlabel('Frequency')
        ylabel('Probabilidy Density')
        box off
        
        %save density estimation
        xdata_hc{par} = r1{1}.XData; ydata_hc{par} = r1{1}.YData;
        xdata_atyp{par} = r2{1}.XData; ydata_atyp{par} = r2{1}.YData;
        xdata_pd{par} = r3{1}.XData; ydata_pd{par} = r3{1}.YData;
        
        %some jitter for groups & Amplitudes
        jit = -.21 : 0.01 : +.2;
        
        jit_hc_amp   = jit(randi([1,numel(jit)],length(Amp_hc),1))';
        jit_atyp_amp = jit(randi([1,numel(jit)],length(Amp_atyp),1))';
        jit_pd_amp   = jit(randi([1,numel(jit)],length(Amp_pd),1))';
        
        h3 = subplot(3,2,6);
        scatter(Amp_hc, 1 + jit_hc_amp,'MarkerEdgeColor',col(1,:),'MarkerFaceColor',col(1,:),'MarkerFaceAlpha',.5);
        hold on
        scatter(Amp_atyp, 2 + jit_atyp_amp,'MarkerEdgeColor',col(4,:),'MarkerFaceColor',col(4,:),'MarkerFaceAlpha',.5);
        hold on
        scatter(Amp_pd, 3 + jit_pd_amp,'MarkerEdgeColor',col(6,:),'MarkerFaceColor',col(6,:),'MarkerFaceAlpha',.5);
        line([mean(Amp_hc),mean(Amp_hc)],[ 1-0.2, 1+0.2],'Color',[0 0 0],'LineWidth',3)
        line([mean(Amp_atyp),mean(Amp_atyp)],[ 2-0.2, 2+0.2],'Color',[0 0 0],'LineWidth',3)
        line([mean(Amp_pd),mean(Amp_pd)],[ 3-0.2, 3+0.2],'Color',[0 0 0],'LineWidth',3)
        ylim([0,4])
        hold on
        
        ylabel('Groups')
        xlabel('Amplitude')
        
        yticklabels(' ')
        
        title('Amplitude of Oscillations')
        
        legend('Hc','Atyp. PD','Idio. PD')
        
        pause(1)
        
        if strcmpi(use_video,'yes')
            %get frame
            F{par} = getframe(gcf);
        end
        
        if par  < length(ycoords)
            cla(h1)
            cla(h2)
            cla(h3)
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%
%%% RUN THIS ... %%% Peak Frequency
%%%%%%%%%%%%%%%%%%%%

load([mpath,'/parcel/time/hc01/hc01_parcel_time.mat']);
parcel2source.pos = parcel_time.brainordinate.pos * 10;     %units: cm -> mm
parcel2source.inside = parcel_time.brainordinate.mask ~= 0; %inside brain positions
parcel2source.mask = parcel_time.brainordinate.mask;        %mask
parcel2source.unit = 'mm';

if strcmpi(use_amp,'no')
    figure('Position',[50 50 1200 400])
    %group parcels together [force if _L ... appear then also _R and vice versa]
    for par = 1 : length(ycoords)
        
        %store parcels for cluster
        clu = [];
        %apply yPos to Parcel Centroids
        for c = 1:length(ParcelCentroid)
            
            yParcel = unique( ParcelCentroid(c,2) );
            
            %distance (cm) to any of the choose y positions
            distances = abs( yParcel - ycoords(par) ) <= y_margin;
            %check if all three coordinates are in the parcel
            check = sum(sum(distances,2),1) >= min( [3,length(yParcel)] );
            
            if any(check)
                %if a parcel fullfills the 'y - criteria' store its index
                clu = cat(1,clu,c);
                %and take its friend on the other hemisphere
                if endsWith( parcel.masklabel{c},'R')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_R','_L')) ) );
                    clu = cat(1,clu,friend_index);
                elseif endsWith( parcel.masklabel{c},'L')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_L','_R')) ) );
                    clu = cat(1,clu,friend_index);
                end
            end
            
        end
        
        %exclude parcels ...
        clu = unique( clu(~ismember(clu,idx_exclude)) );
        %and in the last-1 step also disregard Frontal_Sup. Looks MUCH better then
        if par == length(ycoords)-1; idx_frontal_sup = [find(contains(parcel.masklabel,'Frontal_Sup_L')),find(contains(parcel.masklabel,'Frontal_Sup_R'))]; clu = clu(~ismember(clu,idx_frontal_sup)) ;end
        
        %Parcel Labels
        parcel2source.reg = ismember(parcel2source.mask,clu);
        
        %highlight region 'reg' on mesh -> on mesh
        cfg = [];
        cfg.parameter = 'reg';
        cfg.downsample = 1;
        cfg.method = 'cubic';
        SourceInterpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
        %highlight the positions for the specific region
        highlight.pos = SourceInterpol.pos(logical(SourceInterpol.reg),:);
        
        h1 = subplot(1,2,1);
        %basic background plot
        ft_plot_mesh( mesh, 'facealpha',0.1 ); hold on;
        %add region in another color
        ft_plot_mesh( highlight, 'vertexcolor',[.8 .8 .8] )
        view([90,0])
        
        %current parcels of interest
        parcels_of_interest = parcel.masklabel(clu);
        
        %some variables that are used for plotting & analysing
        Osc_hc = []; Osc_atyp = []; Osc_pd = [];                                % <- Get Center Frequencies of Oscillations of the group (cumulative over all parcels of interest)
        Amp_hc = []; Amp_atyp = []; Amp_pd = [];                                % <- Get absolute Power / Amplitude of the group (cumulative over all parcels of interest)
        Abs_hc = []; Abs_atyp = []; Abs_pd = [];                                % <- Get absolute number of Oscillations of the group (Number of Oscillations per parcel and Group)
        
        for a = 1:length(parcels_of_interest)
            %area
            area  = find(contains(SpecFeat.labels,parcels_of_interest{a}));
            %oscillation data
            Os = cell(length(groups),1);
            %amplitude data
            Am = cell(length(groups),1);
            
            for n = 1:length(groups)
                
                %choose group
                gr = groups{n};
                
                for k = 1:length(subjects)
                    if contains(subjects(k),gr) && ~contains(cbs_info.(subjects{k}).rest.tremor,'yes')
                        
                        %oscillations in 'area' for subject k
                        osc_candidates = SpecFeat.Osc{k}{area};
                        
                        %amplitude to the oscillations in 'area' for subject k
                        amp_candidates = SpecFeat.Amp{k}{area};
                        
                        if length(intervall) == 2
                            %apply intervall restrictions
                            idx = logical( ( osc_candidates > intervall(1) ) .* ( osc_candidates < intervall(2) ) );
                            %prune osc_candidates
                            osc_candidates = osc_candidates(idx);
                            %prune amp_candidates
                            amp_candidates = amp_candidates(idx);
                        end
                        
                        if strcmpi(use_amp,'yes')
                            %select the major ('maximal') oscillation
                            [amp_candidates,ch] = max(amp_candidates);
                            %choose oscillation
                            osc_candidates = osc_candidates(ch);
                        end
                        
                        if isrow(osc_candidates)
                            osc_candidates = osc_candidates';
                        end
                        if isrow(amp_candidates)
                            amp_candidates = amp_candidates';
                        end
                        
                        %save oscillations and amplitudes
                        Os{n} = vertcat(Os{n},osc_candidates);
                        Am{n} = vertcat(Am{n},amp_candidates);
                        
                    end
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Oscillatory Center Frequency %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Absolute number of oscillations in this area and intervall for the different groups
            Abs_hc(a) = length(Os{1});
            Abs_atyp(a) = length(Os{2});
            Abs_pd(a) = length(Os{3});
            
            %Center Frequency of Oscillations
            Osc_hc  = vertcat( Osc_hc,Os{1} );
            Osc_atyp= vertcat( Osc_atyp,Os{2} );
            Osc_pd  = vertcat( Osc_pd,Os{3} );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Amplitudes in defined Intervall %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %get Power for Parcel 'a' for all three groups as three elements of the array Pow
            Pow = cellfun(@sum,Am);
            
            %Amplitude of Center Frequency
            Amp_hc  = vertcat( Amp_hc,Am{1} );
            Amp_atyp= vertcat( Amp_atyp,Am{2} );
            Amp_pd  = vertcat( Amp_pd,Am{3} );
            
            %Average Amplitude per parcel and subject as given by group identification
            Avg_Pow_hc  = sum(Amp_hc) / (a * sum(indices.hc));
            Avg_Pow_atyp= sum(Amp_atyp) / (a * sum(indices.aps));
            Avg_Pow_pd  = sum(Amp_pd) /(a * sum(indices.ips));
            
            %clear Pow for security reasons ... ;)
            clear Pow
            
        end
        
        %Average Number of Oscillations per parcel and subject as given by group identification
        Avg_N_Osc_hc  = sum(Abs_hc)/ (length(parcels_of_interest) * sum(indices.hc));
        Avg_N_Osc_atyp = sum(Abs_atyp)/ (length(parcels_of_interest) * sum(indices.aps));
        Avg_N_Osc_pd  = sum(Abs_pd)/(length(parcels_of_interest) * sum(indices.ips));
        
        %store all Oscillations of Region a
        OscRegionHc{par}= Osc_hc;
        OscRegionAtyp{par}= Osc_atyp;
        OscRegionPd{par}= Osc_pd;
        
        %store Group Median Oscillations (for later Plotting)
        OscMed = vertcat(OscMed,[median(Osc_hc),median(Osc_atyp),median(Osc_pd)]);
        OscIQR = vertcat(OscIQR,[quantile(Osc_hc,.25),quantile(Osc_hc,.75),quantile(Osc_atyp,.25),quantile(Osc_atyp,.75),quantile(Osc_pd,.25),quantile(Osc_pd,.75)]);
        %store Group Mean Oscillations (for later Plotting)
        OscAvg = vertcat(OscAvg,[mean(Osc_hc),mean(Osc_atyp),mean(Osc_pd)]);
        OscStd = vertcat(OscStd,[std(Osc_hc),std(Osc_atyp),std(Osc_pd)]);
        
        %store Group Amplitude (for later Plotting)
        AmpMed = vertcat(AmpMed,[median(Amp_hc),median(Amp_atyp),median(Amp_pd)]);
        AmpAvg = vertcat(AmpAvg,[mean(Amp_hc),mean(Amp_atyp),mean(Amp_pd)]);
        AmpIQR = vertcat(AmpIQR,[quantile(Amp_hc,.25),quantile(Amp_hc,.75),quantile(Amp_atyp,.25),quantile(Amp_atyp,.75),quantile(Amp_pd,.25),quantile(Amp_pd,.75)]);
        %Absolute number of oscillations over y-axis
        Abs_avg_All = vertcat(Abs_avg_All,[Avg_N_Osc_hc,Avg_N_Osc_atyp,Avg_N_Osc_pd]);
        
        %title('N Oscillations')
        
        h2 = subplot(1,2,2);
        r1 = raincloud_plot(Osc_hc, 'box_on', 1, 'color', col(1,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.12, 'dot_dodge_amount', .12,'box_col_match', 0 );
        r2 = raincloud_plot(Osc_atyp, 'box_on', 1, 'color', col(4,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.32, 'dot_dodge_amount', .32,'box_col_match', 0 );
        r3 = raincloud_plot(Osc_pd, 'box_on', 1, 'color', col(6,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.52, 'dot_dodge_amount', .52,'box_col_match', 0 );
        legend([r1{1} r2{1} r3{1}], {'HC', 'APS', 'IPS'})
        title('Center of Frequency')
        set(gca,'XLim', [1,35]);
        set(gca,'YLim', [-0.065,0.1]);
        xlabel('Frequency')
        ylabel('Probabilidy Density')
        box off
        
        %save density estimation
        xdata_hc{par} = r1{1}.XData; ydata_hc{par} = r1{1}.YData;
        xdata_atyp{par} = r2{1}.XData; ydata_atyp{par} = r2{1}.YData;
        xdata_pd{par} = r3{1}.XData; ydata_pd{par} = r3{1}.YData;
        
        pause(1)
        
        if strcmpi(use_video,'yes')
            %get frame
            F{par} = getframe(gcf);
        end
        
        if par  < length(ycoords)
            cla(h1)
            cla(h2)
        end
        
    end
elseif strcmpi(use_amp,'yes')
    figure('Position',[50 50 1200 400])
    %group parcels together [force if _L ... appear then also _R and vice versa]
    for par = 1 : length(ycoords)
        
        %store parcels for cluster
        clu = [];
        %apply yPos to Parcel Centroids
        for c = 1:length(ParcelCentroid)
            
            yParcel = unique( ParcelCentroid(c,2) );
            
            %distance (cm) to any of the choose y positions
            distances = abs( yParcel - ycoords(par) ) <= y_margin;
            %check if all three coordinates are in the parcel
            check = sum(sum(distances,2),1) >= min( [3,length(yParcel)] );
            
            if any(check)
                %if a parcel fullfills the 'y - criteria' store its index
                clu = cat(1,clu,c);
                %and take its friend on the other hemisphere
                if endsWith( parcel.masklabel{c},'R')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_R','_L')) ) );
                    clu = cat(1,clu,friend_index);
                elseif endsWith( parcel.masklabel{c},'L')
                    friend_index = find( strcmpi( parcel.masklabel, (replace(parcel.masklabel{c},'_L','_R')) ) );
                    clu = cat(1,clu,friend_index);
                end
            end
            
        end
        
        %exclude parcels ...
        clu = unique( clu(~ismember(clu,idx_exclude)) );
        %and in the last-1 step also disregard Frontal_Sup. Looks MUCH better then
        if par == length(ycoords)-1; idx_frontal_sup = [find(contains(parcel.masklabel,'Frontal_Sup_L')),find(contains(parcel.masklabel,'Frontal_Sup_R'))]; clu = clu(~ismember(clu,idx_frontal_sup)) ;end
        
        %Parcel Labels
        parcel2source.reg = ismember(parcel2source.mask,clu);
        
        %highlight region 'reg' on mesh -> on mesh
        cfg = [];
        cfg.parameter = 'reg';
        cfg.downsample = 1;
        cfg.method = 'cubic';
        SourceInterpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
        %highlight the positions for the specific region
        highlight.pos = SourceInterpol.pos(logical(SourceInterpol.reg),:);
        
        h1 = subplot(1,2,1);
        %basic background plot
        ft_plot_mesh( mesh, 'facealpha',0.1 ); hold on;
        %add region in another color
        ft_plot_mesh( highlight, 'vertexcolor',[.8 .8 .8] )
        view([90,0])
        
        %current parcels of interest
        parcels_of_interest = parcel.masklabel(clu);
        
        %some variables that are used for plotting & analysing
        Osc_hc = []; Osc_atyp = []; Osc_pd = [];                                % <- Get Center Frequencies of Oscillations of the group (cumulative over all parcels of interest)
        Amp_hc = []; Amp_atyp = []; Amp_pd = [];                                % <- Get absolute Power / Amplitude of the group (cumulative over all parcels of interest)
        Abs_hc = []; Abs_atyp = []; Abs_pd = [];                                % <- Get absolute number of Oscillations of the group (Number of Oscillations per parcel and Group)
        
        for a = 1:length(parcels_of_interest)
            %area
            area  = find(contains(SpecFeat.labels,parcels_of_interest{a}));
            %oscillation data
            Os = cell(length(groups),1);
            %amplitude data
            Am = cell(length(groups),1);
            
            for n = 1:length(groups)
                
                %choose group
                gr = groups{n};
                
                for k = 1:length(subjects)
                    if contains(subjects(k),gr) && ~contains(cbs_info.(subjects{k}).rest.tremor,'yes')
                        
                        %oscillations in 'area' for subject k
                        osc_candidates = SpecFeat.Osc{k}{area};
                        
                        %amplitude to the oscillations in 'area' for subject k
                        amp_candidates = SpecFeat.Amp{k}{area};
                        
                        if length(intervall) == 2
                            %apply intervall restrictions
                            idx = logical( ( osc_candidates > intervall(1) ) .* ( osc_candidates < intervall(2) ) );
                            %prune osc_candidates
                            osc_candidates = osc_candidates(idx);
                            %prune amp_candidates
                            amp_candidates = amp_candidates(idx);
                        end
                        
                        if strcmpi(use_amp,'yes')
                            %select the major ('maximal') oscillation
                            [amp_candidates,ch] = max(amp_candidates);
                            %choose oscillation
                            osc_candidates = osc_candidates(ch);
                        end
                        
                        if isrow(osc_candidates)
                            osc_candidates = osc_candidates';
                        end
                        if isrow(amp_candidates)
                            amp_candidates = amp_candidates';
                        end
                        
                        %save oscillations and amplitudes
                        Os{n} = vertcat(Os{n},osc_candidates);
                        Am{n} = vertcat(Am{n},amp_candidates);
                        
                    end
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Oscillatory Center Frequency %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Absolute number of oscillations in this area and intervall for the different groups
            Abs_hc(a) = length(Os{1});
            Abs_atyp(a) = length(Os{2});
            Abs_pd(a) = length(Os{3});
            
            %Center Frequency of Oscillations
            Osc_hc  = vertcat( Osc_hc,Os{1} );
            Osc_atyp= vertcat( Osc_atyp,Os{2} );
            Osc_pd  = vertcat( Osc_pd,Os{3} );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Amplitudes in defined Intervall %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %get Power for Parcel 'a' for all three groups as three elements of the array Pow
            Pow = cellfun(@sum,Am);
            
            %Amplitude of Center Frequency
            Amp_hc  = vertcat( Amp_hc,Am{1} );
            Amp_atyp= vertcat( Amp_atyp,Am{2} );
            Amp_pd  = vertcat( Amp_pd,Am{3} );
            
            %Average Amplitude per parcel and subject as given by group identification
            Avg_Pow_hc  = sum(Amp_hc) / (a * sum(indices.hc));
            Avg_Pow_atyp= sum(Amp_atyp) / (a * sum(indices.aps));
            Avg_Pow_pd  = sum(Amp_pd) /(a * sum(indices.ips));
            
            %clear Pow for security reasons ... ;)
            clear Pow
            
        end
                
        %store all Oscillations of Region a
        OscRegionHc{par}= Osc_hc;
        OscRegionAtyp{par}= Osc_atyp;
        OscRegionPd{par}= Osc_pd;
        
        %store Group Median Oscillations (for later Plotting)
        OscMed = vertcat(OscMed,[median(Osc_hc),median(Osc_atyp),median(Osc_pd)]);
        OscIQR = vertcat(OscIQR,[quantile(Osc_hc,.25),quantile(Osc_hc,.75),quantile(Osc_atyp,.25),quantile(Osc_atyp,.75),quantile(Osc_pd,.25),quantile(Osc_pd,.75)]);
        %store Group Mean Oscillations (for later Plotting)
        OscAvg = vertcat(OscAvg,[mean(Osc_hc),mean(Osc_atyp),mean(Osc_pd)]);
        OscStd = vertcat(OscStd,[std(Osc_hc),std(Osc_atyp),std(Osc_pd)]);
        
        %store Group Amplitude (for later Plotting)
        AmpMed = vertcat(AmpMed,[median(Amp_hc),median(Amp_atyp),median(Amp_pd)]);
        AmpAvg = vertcat(AmpAvg,[mean(Amp_hc),mean(Amp_atyp),mean(Amp_pd)]);
        AmpIQR = vertcat(AmpIQR,[quantile(Amp_hc,.25),quantile(Amp_hc,.75),quantile(Amp_atyp,.25),quantile(Amp_atyp,.75),quantile(Amp_pd,.25),quantile(Amp_pd,.75)]);
        
        h2 = subplot(1,2,2);
        r1 = raincloud_plot(Osc_hc, 'box_on', 1, 'color', col(1,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.12, 'dot_dodge_amount', .12,'box_col_match', 0 );
        r2 = raincloud_plot(Osc_atyp, 'box_on', 1, 'color', col(4,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.32, 'dot_dodge_amount', .32,'box_col_match', 0 );
        r3 = raincloud_plot(Osc_pd, 'box_on', 1, 'color', col(6,:), 'alpha', 0.4,...
            'box_dodge', 1, 'box_dodge_amount', 0.52, 'dot_dodge_amount', .52,'box_col_match', 0 );
        legend([r1{1} r2{1} r3{1}], {'HC', 'APS', 'IPS'})
        title('Center of Frequency')
        set(gca,'XLim', [1,35]);
        %set(gca,'YLim', [-0.07,0.1]);
        xlabel('Frequency')
        ylabel('Probabilidy Density')
        box off
        
        %save density estimation
        xdata_hc{par} = r1{1}.XData; ydata_hc{par} = r1{1}.YData;
        xdata_atyp{par} = r2{1}.XData; ydata_atyp{par} = r2{1}.YData;
        xdata_pd{par} = r3{1}.XData; ydata_pd{par} = r3{1}.YData;
                
        pause(1)
        
        if strcmpi(use_video,'yes')
            %get frame
            F{par} = getframe(gcf);
        end
        
        if par  < length(ycoords)
            cla(h1)
            cla(h2)
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%
%%% END OF VIDEO %%%
%%%%%%%%%%%%%%%%%%%%

if strcmpi(use_video,'yes')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%___ Create Movie ___%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %make directory
    if ~exist([mpath,'/parcel/ana_shift/brainmovie'],'dir')
        mkdir([mpath,'/parcel/ana_shift/brainmovie'])
    end
    
    %initiate VideoObject
    writerObj = VideoWriter([mpath,'/parcel/ana_shift/brainmovie/SpecFeat_GroupVideo_UseAmp_',use_amp,'.avi']);
    writerObj.FrameRate = 1;
    %open the video writer
    open(writerObj);
    %write the frames to the video
    for l=1:length(F)
        %convert the image to a frame
        frame = F{l} ;
        writeVideo(writerObj, frame);
    end
    %close the writer object
    close(writerObj);
    
end

%Summary Figure
figure('Position',[50 50 750 600])
%figure('Position',[25 25 700 300])
subplot(2,2,1)
plot(ycoords,OscMed(:,1),'Color',col(1,:),'LineWidth',3)
hold on
plot(ycoords,OscMed(:,2),'Color',col(4,:),'LineWidth',3)
hold on
plot(ycoords,OscMed(:,3),'Color',col(6,:),'LineWidth',3)
xlim([min(ycoords),max(ycoords)])
ylim([ceil(min(min(OscIQR))-2),ceil(max(max(OscIQR))+2)])
title('Center Frequency [Median]')
xlabel('y-coordinates')
ylabel('Center Frequency')
legend({'Hc','Atyp. PD','Idio. PD'},'FontSize',11,'FontWeight','bold')

subplot(2,2,2)
plot(ycoords,AmpMed(:,1),'Color',col(1,:),'LineWidth',3)
hold on
plot(ycoords,AmpMed(:,2),'Color',col(4,:),'LineWidth',3)
hold on
plot(ycoords,AmpMed(:,3),'Color',col(6,:),'LineWidth',3)
xlim([min(ycoords),max(ycoords)])
ylim([0,max(max(AmpAvg))+0.1])
title('Amplitude [Median]')
xlabel('y-coordinates')
ylabel('Amplitude')
legend({'Hc','Atyp. PD','Idio. PD'},'FontSize',11,'FontWeight','bold')

saveas(gcf,[mpath,'/parcel/ana_shift/Oscillations_summary_AmpUsed_',use_amp,'.png'])

close

%%%%%%%%%%%%%%%%%%%%%%
%%% Save Variables %%%
%%%%%%%%%%%%%%%%%%%%%%

%create group structure with content for statistical testing
GroupRegionOsc.OscRegionHc = OscRegionHc;
GroupRegionOsc.OscRegionPd = OscRegionPd;
GroupRegionOsc.OscRegionAtyp = OscRegionAtyp;
%y-coordinates
GroupRegionOsc.ycoords = ycoords;

%save structure
save([mpath,'/parcel/ana_shift/GroupRegionOsc.mat'],'GroupRegionOsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Density over ycoords %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(use_amp,'no')
    
    figure('Position',[-125,50,1500,500])
    col_hc = cbrewer('seq','Greens',length(ycoords) + 5); col_hc = col_hc(5 + 1:length(ycoords) + 5,:);
    col_atyp = cbrewer('seq','Reds',length(ycoords) + 5); col_atyp = col_atyp(5 + 1:length(ycoords) + 5,:);
    col_pd = cbrewer('seq','Oranges',length(ycoords) + 5); col_pd = col_pd(5 + 1:length(ycoords) + 5,:);
    %first subplot
    for k = 1:length(ycoords)
        subplot(1,3,1)
        plot3(xdata_hc{k},repmat(ycoords(k),length(xdata_hc{k}),1),ydata_hc{k},'Color',col_hc(k,:),'LineWidth',4)
        hold on
    end
    %change settings
    title('HC','FontSize',22); xlim([0 35]); view([20 40])
    %changes axis
    ax = gca;
    ax.FontSize = 20;
    %change y axis
    ax.YLabel.String = 'y-axis';
    ax.YLabel.Rotation = 66;
    ax.YLabel.FontSize = 22;
    ax.YLabel.Position = [38, -1.8, -0.02];
    ax.YLim = [-10,7];
    %change x axis
    ax.XLabel.String = 'Frequency';
    ax.XLabel.Rotation = -10;
    ax.XLabel.FontSize = 22;
    ax.XLabel.Position = [22, -10, -0.02];
    %change z axis
    ax.ZLabel.String = 'Density';
    ax.ZLabel.Rotation = 0;
    ax.ZLabel.FontSize = 22;
    ax.ZLabel.Position = [-21,13,0];
    ax.ZLim = [0 0.115];
    %second subplot
    for k = 1:length(ycoords)
        subplot(1,3,2)
        plot3(xdata_pd{k},repmat(ycoords(k),length(xdata_pd{k}),1),ydata_pd{k},'Color',col_pd(k,:),'LineWidth',4)
        hold on
    end
    title('PD','FontSize',16); xlim([0 35]); view([20 40])
    %changes axis
    ax = gca;
    ax.FontSize = 20;
    %change y axis
    ax.YLabel.String = 'y-axis';
    ax.YLabel.Rotation = 66;
    ax.YLabel.FontSize = 22;
    ax.YLabel.Position = [38, -1.8, -0.02];
    ax.YLim = [-10,7];
    %change x axis
    ax.XLabel.String = 'Frequency';
    ax.XLabel.Rotation = -10;
    ax.XLabel.FontSize = 22;
    ax.XLabel.Position = [22, -10, -0.02];
    %change z axis
    ax.ZLabel.String = 'Density';
    ax.ZLabel.Rotation = 0;
    ax.ZLabel.FontSize = 22;
    ax.ZLabel.Position = [-21,13,0];
    ax.ZLim = [0 0.115];
    %third subplot
    for k = 1:length(ycoords)
        subplot(1,3,3)
        plot3(xdata_atyp{k},repmat(ycoords(k),length(xdata_atyp{k}),1),ydata_atyp{k},'Color',col_atyp(k,:),'LineWidth',4)
        hold on
    end
    title('APS');
    xlim([0 35]);
    view([20 40])
    %changes axis
    ax = gca;
    ax.FontSize = 20;
    %change y axis
    ax.YLabel.String = 'y-axis';
    ax.YLabel.Rotation = 66;
    ax.YLabel.FontSize = 22;
    ax.YLabel.Position = [38, -1.8, -0.02];
    ax.YLim = [-10,7];
    %change x axis
    ax.XLabel.String = 'Frequency';
    ax.XLabel.Rotation = -10;
    ax.XLabel.FontSize = 22;
    ax.XLabel.Position = [22, -10, -0.02];
    %change z axis
    ax.ZLabel.String = 'Density';
    ax.ZLabel.Rotation = 0;
    ax.ZLabel.FontSize = 22;
    ax.ZLabel.Position = [-21,13,0];
    ax.ZLim = [0 0.115];
    
    print([mpath,'/parcel/ana_shift/Oscillations_Densities_AmpUse_',use_amp,'.tiff'],'-dtiff','-r300');
    
    close
    
end