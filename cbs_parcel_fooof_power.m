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

%subjects
subjects = fieldnames(cbs_info);
%catch problematic subjects
catch_sub = [];

%load source information (for positions, unit, labels)
load([mpath,'/parcel/time/hc01/hc01_parcel_time.mat']);

%flip Power & Oscillations (according to UPRDS / Neuropsychological scores) -> right brain hemisphere the more affected one
flipIt = 'yes';

%%%%%%%%%%%%%%%%%%
%%%___script___%%%
%%%%%%%%%%%%%%%%%%
                               
%run the subsequent two lines only if you want to delete APS > 70 & 4 largest UPDRS of CBS
subjects = cbs_clean_subjects(subjects,cbs_info,'exclude',...
                              {'psp03','psp04','psp07','psp12','psp13','psp17','cbs03','cbs06','cbs14','cbs16','cbs18','cbs09','cbs10'});

%tremor subjects (from pd)
tremor = cbs_clean_subjects(subjects,cbs_info,{'rest','tremor'},'yes');
%notremor subjects
notremor = cbs_clean_subjects(subjects,cbs_info,{'rest','tremor'},'-yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%___ Neuropsychological Scores ___%%% -> [more affected bodyside, relative severity, absolute scores]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to xlsx-sheet
testpath = 'C:/Users/mkroe/Documents/Ausbildung/Promotion/project/RestingState/Excel_Subject_Informations/resting_state_study_I/';

%tests to use [which are in the xlsx-file on seperate sheets]
use_test = {'overall'}; %{'updrs_pre_off','tulia','goldenberg','moca'};
use_field = {'laterality'};

%excel structure
excel = struct();

for k = 1:length(use_test)
    excel.(use_test{k}) = readtable([testpath,'COHORT_Tests.xlsx'],'Sheet',use_test{k});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ___Subject Loop___ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%n-trials
ntrials = zeros(1,length(subjects));

%parcel labels
labels = parcel_time.label;

%prepare oscillations (center frequency) and amplitude (power) arrays
osc = cell(1,length(subjects)); for i = 1:length(subjects); osc{i} = cell(1,length(labels)); end
amp = cell(1,length(subjects)); for i = 1:length(subjects); amp{i} = cell(1,length(labels)); end

for i = 1:length(subjects)
    
    try
        
        clear parcel_time
        
        %load subject fooof data
        load([mpath,'/fooof_flexible/',subjects{i},'/results_fooof/parcel/',subjects{i},'_parcel_fooof.mat']);
        
        %"preallocate" fooof_spec. (note it will only be called once. Once the variable fooof_spec exist it the following line will not be executed)
        freq_range = [2,48];
        freq_res   = 0.5;
        
        if ~exist('fooof_spec','var')
            fooof_spec = zeros( length(labels), length( freq_range(1):freq_res:freq_range(2) ), length(subjects) );
        end
        
        %extract fooofed spectrum
        for k = 1:length(labels)
            
            %go on here ....
            fooof_spec(k,:,i) = cbs_extract_fooof_flexible( Parcel.(labels{k}), 'spectrum_flat', freq_range, freq_res );
            
            % get oscillatory peaks and amplitude
            [o,a] = cbs_extract_fooof_oscillations( Parcel.(labels{k}) );
            osc{i}{k} = o(:);
            amp{i}{k} = a(:);
            
        end
        
        %load power spectrum -> only to get n trials
        load([mpath,'/parcel/time/',subjects{i},'/',subjects{i},'_parcel_time.mat'])
        %this line actually is shitty. I should get the n-trials in cbs_cleaning and save it in the .mat. But this I do with the LinuxRechner
        if isfield(parcel_time,'trial')
            ntrials(i) = length( parcel_time.trial );
        end
        
    catch
        catch_sub = horzcat(catch_sub,i);
    end
    
end
warning(['subjects excluded from source analysis: ' strjoin( subjects(catch_sub),' ')])

%index subjects with parcel data (remove the rest)
idx = squeeze( sum(sum(abs(fooof_spec),1),2) ) ~= 0;

%power related data
freqs = freq_range(1) : freq_res : freq_range(2);   %freqs
subjects = subjects(idx);                           %remaining subjects
fooof_spec = fooof_spec(:,:,idx);                   %remaining data

%spectral details data
osc = osc(idx); %clean 'osc' as well (so it matches to the cleaned subjects variable)
amp = amp(idx); %clean 'amp' as well (so it matches to the cleaned subjects variable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flip Osc/Amp & Power %%% <- can be specified at the script beginning. If flipIt = 'yes' the resulting data matrices are used for all further scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%___Oscillations & Amplitude
if strcmpi(flipIt,'yes')
    
    flipped_osc = osc; %structure to reorganize to: right brain hemisphere more affected & left brain hemisphere less affected
    flipped_amp = amp; %structure to reorganize to: right brain hemisphere more affected & left brain hemisphere less affected
    
    %check for subject if activity flipping is necessary
    for s = 1:length(subjects)
        
        %index in excel table
        idx_excel = strcmpi(excel.(use_test{:}).subject_id,subjects(s));
        %index in osc / amp
        idx_data = strcmpi(subjects,subjects(s));
        
        %which body hemisphere is more affected by the disease
        more_affected_bodyside = excel.(use_test{:}).(use_field{:})(idx_excel);
        
        %right brain hemisphere should be the more affected one for all subjects -> so only flip activity if right body hemisphere is affected
        if strcmpi(more_affected_bodyside,'right')
            flipped_osc{idx_data} = cbs_mirror_activity_parcel(labels,osc{idx_data}')';
            flipped_amp{idx_data} = cbs_mirror_activity_parcel(labels,amp{idx_data}')';
        end
        
    end    
    
    osc = flipped_osc;
    amp = flipped_amp;
    
    clear flipped_osc flipped_amp
    
end
%___Power
if strcmpi(flipIt,'yes')
    
    flipped_fooof_spec = fooof_spec; %structure to reorganize to: right brain hemisphere more affected & left brain hemisphere less affected
    
    %check for subject if activity flipping is necessary
    for s = 1:size(subjects,1)
        
        %index in excel table
        idx_excel = strcmpi(excel.(use_test{:}).subject_id,subjects(s));
        %index in fooof_spec
        idx_fooof_spec = strcmpi(subjects,subjects(s));
        
        %which body hemisphere is more affected by the disease
        more_affected_bodyside = excel.(use_test{:}).(use_field{:})(idx_excel);
        
        %right brain hemisphere should be the more affected one for all subjects -> so only flip activity if right body hemisphere is affected
        if strcmpi(more_affected_bodyside,'right')
            flipped_fooof_spec(:,:,idx_fooof_spec) = cbs_mirror_activity_parcel(labels,fooof_spec(:,:,idx_fooof_spec));
        end
        
    end    
    
    fooof_spec = flipped_fooof_spec;
    
    clear flipped_fooof_spec
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Fooof Parcel Structure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%folder to put results
if ~exist([mpath,'/parcel/ana_power'],'dir')
    mkdir([mpath,'/parcel/ana_power'])
end

%Make a struct to save all subjects with parcel data
Fooof.freqs = freqs;
Fooof.sub = subjects;
Fooof.labels = labels;
Fooof.fooof_spec = fooof_spec;
Fooof.ntrials = ntrials';
Fooof.flipped_affected2rightHemisphere = [flipIt,'_',use_test{:},'_',use_field{:}];

if any(4 == freqs(~isnan( sum(sum(fooof_spec,1),3) ))) && any(30 == freqs(~isnan( sum(sum(fooof_spec,1),3) ))) && ~any( isnan(freqs(~isnan( sum(sum(fooof_spec,1),3) ))) )
    %all is good...
    disp('Good. All frequencies needed for further analysis are fine')

    %prune according to frequencies
    Fooof.fooof_spec = Fooof.fooof_spec( :, logical( (Fooof.freqs >= 4) .* (Fooof.freqs <= 30)), : );    
    fooof_spec = fooof_spec( :, logical( (Fooof.freqs >= 4) .* (Fooof.freqs <= 30)), : );
    
    Fooof.freqs = Fooof.freqs(logical( (Fooof.freqs >= 4) .* (Fooof.freqs <= 30) ));
    freqs = freqs( logical( (freqs >= 4) .* (freqs <= 30) ) );
    
else
    error('Something is wrong in the "fooof_spec" variable.')
end

%Save the spectra for inspection in another folder
save([mpath,'/parcel/ana_power/fooof_spectra.mat'],'Fooof');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Spectral Details Structure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%folder to put results
if ~exist([mpath,'/parcel/ana_shift'],'dir')
    mkdir([mpath,'/parcel/ana_shift'])
end

SpecFeat.Osc = osc;
SpecFeat.Amp = amp;
SpecFeat.labels = labels;
SpecFeat.subjects = subjects;
SpecFeat.flipped_affected2rightHemisphere = [flipIt,'_',use_test{:},'_',use_field{:}];

%Save the spectral details for inspection in another folder
save([mpath,'/parcel/ana_shift/spectral_features.mat'],'SpecFeat');

clear Parcel

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ___Parcel Power___ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

parcel2source.pos = parcel_time.brainordinate.pos * 10;     %units: cm -> mm
parcel2source.inside = parcel_time.brainordinate.mask ~= 0; %inside brain positions
parcel2source.mask = parcel_time.brainordinate.mask;        %mask
parcel2source.unit = 'mm';

%mesh in mm
mesh = ft_read_headshape([ft_path,'/template/anatomy/surface_white_both.mat']);

%subject groups structure
indices.hc = contains(subjects,'hc');
indices.cbs = contains(subjects,'cbs');
indices.psp = contains(subjects,'psp');
indices.aps = logical( contains(subjects,'cbs') + contains(subjects,'psp') );
indices.pd = logical( contains(subjects,'pd') .* ismember(subjects,notremor) );

%colormap [for group difference]
Ccomp = cbrewer('div','RdBu',64);
Csing = cbrewer('seq','OrRd',64);

%groupnames
groupnames = fieldnames(indices);

%code group comparison [row vs column]
disp( strcat(num2str((1:length(groupnames))'),': ',groupnames) )
%use indices to groups and compare means
CodeFooof = {1:length(groupnames),[1 2],[1 3],[1 4],[5 4],[5 2],[5 3]};
%frequencies
freqband = {[4 7.5],[8 12.5],[13 19.5], [20 30]};
freqband_name = {'\theta','\alpha','Low {\beta}','High {\beta}'};

%group difference power [all frequencies] for Brain Video
GDPow = [];

for j = 2%1:length(CodeFooof)
    %group code
    code = CodeFooof{j};
    if length(code) ~= 2 %if there are not two groups to compare, then plot each group individually
        for m = 1:length(code)
            if sum( indices.(groupnames{code(m)}) ) == 0
                warning(['No subjects for group: ',groupnames{code(m)}])
                continue
            end
            %group index
            gr = indices.( groupnames{code(m)} );
            %group powspctrm (change later)
            group_powspctrm = fooof_spec;
            %group ntrails
            gr_ntrials = ntrials(gr);
            %extract spectr
            group_powspctrm = squeeze( mean(group_powspctrm(:,:,gr),3,'omitnan') );
            
            %make surface plots
            for f = 1:length(freqband)
                %frequency indices
                freqidx = ismember(freqs,freqband{f}(1):freq_res:freqband{f}(2));
                %frequency band data
                fb_pow = mean(group_powspctrm(:,freqidx),2);
                
                %put fb_pow from parcel to sources
                parcel2source.pow = parcel_time.brainordinate.mask;
                for n = 1:length(fb_pow)
                    %exchange the mask index with the related frequency band power of that mask
                    parcel2source.pow(parcel2source.pow == n) = fb_pow(n);
                end

                %frequency band power -> on mesh
                cfg = [];
                cfg.parameter = 'pow';
                cfg.downsample = 2;
                cfg.method = 'cubic';
                SourceInterpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
                %make plot
                cfg = [];
                cfg.method = 'surface';
                cfg.funparameter = 'pow';
                cfg.projmethod = 'nearest';
                cfg.camlight = 'no';
                cfg.funcolormap = Csing;
                ft_sourceplot(cfg,SourceInterpol)
                gr_name = replace(groupnames{code(m)},'_','. ');
                c = colorbar;
                c.FontSize = 18;
                c.FontWeight = 'bold';
                c.Label.String = 'log10(Power)';
                c.Label.Position = [-1.5 mean(c.Limits) 0];
                view([90 0])
                
                %save image (right view)
                saveas(gcf,[mpath,'/parcel/ana_power/',groupnames{code(m)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_right.tiff']);                
                %change views and make images (right -> left -> top)
                view([-90 0])
                colorbar off
                saveas(gcf,[mpath,'/parcel/ana_power/',groupnames{code(m)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_left.png']);
                %other side ...
                view([0 90])
                current_fb = freqband_name{f};
                title([current_fb,': ',replace([upper(groupnames{code(m)}(1)),groupnames{code(m)}(2:end)],'_',' ')],'FontSize',18);    
                saveas(gcf,[mpath,'/parcel/ana_power/',groupnames{code(m)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_top.png']);

                close
            end
        end
        
    else %otherwise we have two groups
        %group index
        gr1 = indices.( groupnames{code(1)} );
        gr2 = indices.( groupnames{code(2)} );
        %group ntrails
        gr1_ntrials = ntrials(gr1);
        gr2_ntrials = ntrials(gr2);
        %group1 spectrum: calculate mean and make use of the merging matrix
        group1_powspctrm = squeeze( mean(fooof_spec(:,:,gr1),3,'omitnan') );
        %group2 spectrum: calculate mean and make use of the merging matrix
        group2_powspctrm = squeeze( mean(fooof_spec(:,:,gr2),3,'omitnan') );
        %group difference (all frequencies)
        groupdiff_powspctrm = group1_powspctrm - group2_powspctrm;
        %save group difference in GDPow
        GDPow = cat(3,GDPow,groupdiff_powspctrm);
        
        %make surface plots
        for f = 1:length(freqband)
            %frequency indices
            freqidx = ismember(freqs,freqband{f}(1):freq_res:freqband{f}(2));
            %frequency band data
            fb_pow = mean(groupdiff_powspctrm(:,freqidx),2);
            
            %put fb_pow from parcel to sources
            parcel2source.pow = parcel_time.brainordinate.mask;
            for n = 1:length(fb_pow)
                %exchange the mask index with the related frequency band power of that mask
                parcel2source.pow(parcel2source.pow == n) = fb_pow(n);
            end
            
            %frequency band power -> on mesh
            cfg = [];
            cfg.parameter = 'pow';
            cfg.downsample = 2;
            cfg.method = 'cubic';
            SourceInterpol = ft_sourceinterpolate(cfg,parcel2source,mesh);
            %make plot
            cfg = [];
            cfg.method = 'surface';
            cfg.funparameter = 'pow';
            cfg.projmethod = 'nearest';
            cfg.camlight = 'no';
            cfg.funcolormap = Ccomp;
            cfg.funcolorlim = [-0.2,0.2];
            cfg.colorbar = 'no';
            ft_sourceplot(cfg,SourceInterpol)
            gr1_name = replace(groupnames{code(1)},'_','. ');
            gr2_name = replace(groupnames{code(2)},'_','. ');
            set(gcf,'Position',[10 10 650 700])
            %right view
            view([90 0])
            current_fb = freqband_name{f};
            title([current_fb,': ',upper(gr1_name) ,' vs. ',upper(gr2_name)],'FontSize',32);
            c = colorbar('SouthOutside');
            c.FontSize = 32;
            c.FontWeight = 'bold';
            c.Label.String = {[upper(gr1_name),' < ',upper(gr2_name),'                   ',upper(gr2_name),' < ',upper(gr1_name)],...
                [''],...
                [''],...
                ['log10(Power)']};
            c.Label.Position = [0 4 0];
            print([mpath,'/parcel/ana_power/',groupnames{code(1)},'VS',groupnames{code(2)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_right.tiff'],'-dtiff','-r300');
            %left view
            view([-90 0])
            print([mpath,'/parcel/ana_power/',groupnames{code(1)},'VS',groupnames{code(2)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_left.tiff'],'-dtiff','-r300');
            %change views and make images (right -> left -> top)
            view([0 90])
            print([mpath,'/parcel/ana_power/',groupnames{code(1)},'VS',groupnames{code(2)},'_',num2str(freqband{f}(1)),'to',num2str(freqband{f}(2)),'hz_top.tiff'],'-dtiff','-r300');

            close
            
        end
    end
end
