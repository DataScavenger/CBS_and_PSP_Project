%% Statistical Testing on Parcels of Frequency Bands [Cluster Permutation Test]

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

%load frequency band data
load([mpath,'/parcel/ana_power/fooof_spectra.mat'],'Fooof');

%notremor subjects
notremor = cbs_clean_subjects(Fooof.sub,cbs_info,{'rest','tremor'},'-yes');

%subject groups structure
indices.hc = contains(Fooof.sub,'hc');
indices.cbs = contains(Fooof.sub,'cbs');
indices.psp = contains(Fooof.sub,'psp');
indices.aps = logical( contains(Fooof.sub,'cbs') + contains(Fooof.sub,'psp') );
indices.pd = logical( contains(Fooof.sub,'pd') .* ismember(Fooof.sub,notremor) );

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
groupnames = fieldnames(indices);

%code group comparison (statistics)
disp( strcat(num2str((1:length(groupnames))'),': ',groupnames) )
%use indices to groups and compare means
CodeA = {[1 3],[1 2],[1 4],[1 5],[2 5],[3 5],[4 5]};

%frequencies
freqband = {[4 7.5],[8 12.5],[13 19.5],[20 30]};
freqband_name = {'\theta','\alpha','Low {\beta}','High {\beta}','{\beta}'};

%color for statistical value
Cstat = cbrewer('div','PRGn',64);

%save p-values
p_v = nan(length(CodeA),length(freqband));

%group comparison
for f = 1:length(freqband)
    
    for j = 1:length(CodeA)
        
        %group code
        code = CodeA{j};
        %group index
        gr1 = indices.( groupnames{code(1)} );
        gr2 = indices.( groupnames{code(2)} );
        %freqband from -> to
        from = find(ismember( Fooof.freqs, freqband{f}(1) ));
        to = find(ismember( Fooof.freqs, freqband{f}(2) ));
        %group Frequency Band
        gr1_freq = squeeze( mean( Fooof.fooof_spec(:,from:to,gr1),2 ) );
        gr2_freq = squeeze( mean( Fooof.fooof_spec(:,from:to,gr2),2 ) );
        
        %labels
        labels = Fooof.labels;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% Statistical Test %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        %hack into a fieldtrip-like structure
        gr1_data = {};
        tmp = [];
        for k = 1:size(gr1_freq,2)
            tmp.Pow = squeeze(gr1_freq(:,k));
            tmp.label = labels;
            %act as if it was 'power'
            tmp.freq = 5;
            tmp.dimord = 'chan_freq';
            %put tmp in data
            gr1_data{k} = tmp;
        end
        %hack into a fieldtrip-like structure
        gr2_data = {};
        tmp = [];
        for k = 1:size(gr2_freq,2)
            tmp.Pow = squeeze(gr2_freq(:,k));
            tmp.label = labels;
            %act as if it was 'power'
            tmp.freq = 5;
            tmp.dimord = 'chan_freq';
            %put tmp in data
            gr2_data{k} = tmp;
        end
        
        %___ Statistical Test [Frequency Band Power]
        
        %group frequencies (for statistical testing)
        cfg = [];
        cfg.keepindividual = 'yes';
        cfg.parameter = 'Pow';
        gr1_data = ft_freqgrandaverage(cfg,gr1_data{:});
        gr2_data = ft_freqgrandaverage(cfg,gr2_data{:});
        
        %cluster permutation test
        cfg = [];
        cfg.neighbours = neigh;
        cfg.minnbchan = 1;
        cfg.method = 'montecarlo';
        cfg.statistic = 'ft_statfun_indepsamplesT';
        cfg.parameter = 'Pow';
        cfg.correctm = 'cluster';
        cfg.numrandomization = 10000;
        cfg.alpha = 0.025;
        cfg.clusteralpha = 0.025;
        cfg.clusterstatistic = 'maxsum';
        cfg.tail = 0;
        %statistical parameter settings
        design(1,:) = [ones(1,size(gr1_data.Pow,1))*1,ones(1,size(gr2_data.Pow,1))*2];
        cfg.design = design;
        cfg.ivar = 1;
        %perform statistical test
        s = ft_freqstatistics(cfg,gr1_data,gr2_data);
        %clear 'design'
        clear design
                
        if sum(s.mask) ~= 0
            
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
            
            %save in folder
            if ~isdir([mpath,'/parcel/ana_power/statistic'])
                mkdir([mpath,'/parcel/ana_power/statistic'])
            end
            
            %plot t-values
            cfg = [];
            cfg.method = 'surface';
            cfg.funparameter = 'tval';
            cfg.projmethod = 'nearest';
            cfg.funcolormap = Cstat;
            cfg.camlight = 'no';
            cfg.funcolorlim = [-3,3];
            ft_sourceplot(cfg,tmp_interpol)
            set(gcf,'Position',[10 10 650 700])
            title([freqband_name{f},': ',upper(replace((groupnames{code(1)}),'_','. ')),' vs ',upper(replace((groupnames{code(2)}),'_','. '))],'FontSize',32);
            colorbar off
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
            print([mpath,'/parcel/ana_power/statistic/',groupnames{code(1)},'VS',groupnames{code(2)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_left.tiff'],'-dtiff','-r300');
            %right view
            view([90 0])
            print([mpath,'/parcel/ana_power/statistic/',groupnames{code(1)},'VS',groupnames{code(2)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_right.tiff'],'-dtiff','-r300');
            %change views and make images (right -> left -> top)
            view([0 90])
            print([mpath,'/parcel/ana_power/statistic/',groupnames{code(1)},'VS',groupnames{code(2)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_top.tiff'],'-dtiff','-r300');

            close
            
            p_v(j,f) = min(s.prob);
            
        else
            warning(['No significant clusters | ',replace(groupnames{code(1)},'_',' '),' vs ',replace(groupnames{code(2)},'_',' '),' | ',num2str(freqband{f}(1)),' to ',num2str(freqband{f}(2)),' Hz'])
                        
        end
        
    end
    
end