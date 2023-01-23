%% Time & Frequency Domain Source Reconstruction [LCMV Beamformer]

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

%'catch' subjects
catch_sub = [];

for i = 1:length(subjects)
    
    clearvars hdm data_select_clean spatial_filt trial_source CovData parcel_time parcel_time_seg
    
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ___Load Variables___ %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~exist([mpath,'/headmodels/',subjects{i},'/',subjects{i},'_forward_model.mat'],'file')
            continue
        end
        
        %load headmodel & grid (subject specific)
        load([mpath,'/headmodels/',subjects{i},'/',subjects{i},'_forward_model.mat']);
        
        %load the subject data (time domain)
        load([mpath,'/clean_data/',subjects{i},'_rest_time.mat']);
        
        %the field 'elec' was artificially added in cbs_prep_cleaning by the function ft_rejectcomponent() even though we only select megplanars. It wrongly put the assumption that all MEG channels are eeg channels, which is wrong. The other datasets instead do no contain this. Thus, here it is removed
        if strcmpi(subjects{i},'pd59'); data_select_clean = rmfield(data_select_clean,'elec'); end
        if strcmpi(subjects{i},'hc11'); data_select_clean = rmfield(data_select_clean,'elec'); end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%___ Source Time Computation ___%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        cfg = [];
        cfg.length = 2;
        %cfg.overlap = 0.5;
        CovData = ft_redefinetrial(cfg,data_select_clean);
        
        %change distorted time-field to calculate the covariance matrix (segmented data is to large to be processed)
        CovData.time(:) = CovData.time(1);
        
        %we need the covariance matrix for the calculation of spatial filters (for source reconstruction, which is slightly different from cov(data) for whatever reason)
        cfg = [];
        cfg.covariance = 'yes';
        avg = ft_timelockanalysis(cfg,CovData);
        
        %compute spatial filter
        cfg = [];
        cfg.method = 'lcmv';
        cfg.lcmv.lambda = '5%';
        cfg.headmodel = hdm;
        cfg.sourcemodel = grid;
        cfg.normalize = 'yes'; %normalize the leadfield (scales the spatial filter at a grid point)
        if strcmpi(subjects{i},'cbs01'); load([mpath,'/new_grad_cbs01.mat']); cfg.grad = cbs01_grad_new; end
        if strcmpi(subjects{i},'pd59'); load([mpath,'/new_grad_pd59.mat']); cfg.grad = pd59_grad_new; end        
        cfg.lcmv.projectmom = 'yes';
        cfg.lcmv.reducerank = 2;
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.projectnoise = 'yes';
        source_time_tmp = ft_sourceanalysis(cfg,avg);
        %computation of virtual channels (svd happens implicitely with projectmom = 'yes'. Look up under ft_inverse_lcmv which is run in ft_sourceanalysis)
        
        %Adapt size(filter) to match .pos [These must match, that the script does not crash. Yet, these two lines of code are not important for the used datasets]
        source_time_tmp.avg.filter(cellfun('isempty',source_time_tmp.avg.filter)) = {zeros(1,length(source_time_tmp.avg.label))};
        spatial_filt = cell2mat(source_time_tmp.avg.filter); %spatial filter for subject{i}, Sources X Channels
        
        %Security check. Zero-line sources are not part of a parcel which is needed later
        if any( sum( abs( spatial_filt( parcel.mask(:) > 0 ,: ) ) ,2 ) == 0 )
            error(['At least one source that contributes to the calculation of parcel activation contain zero-line channel-weights,' newline 'which results in a zero-line time-course. However this source should not contribute to the parcel.'])
        end
        
        %Prepare Filter application
        trial_source = cell(1,length(data_select_clean.trial));
        
        %get the data [with all trials] projected into source space with the individual spatial filter
        for j = 1:length(data_select_clean.trial)
            trial_source{j} = spatial_filt * data_select_clean.trial{j}; %Sources X Timepoint X Trails
        end
        
        %we need labels for later analysis (like ft_freqanalysis())
        label = cell(size(trial_source,1),1);
        for s = 1:size(trial_source{1},1)
            label{s} = horzcat('source_',num2str(s));
        end
        
        %save results in a structure
        source_time.pos = parcel.pos;                   %go back to the template_grid and its .pos. source_time_tmp.pos would be the the subject-specific grid which is not appropriate for group comparison
        source_time.trial = trial_source;               %#sources X [sampling_rate*trials]
        source_time.label = label;                      %source label
        source_time.time = data_select_clean.time;      %time points
        source_time.unit = 'cm';
        
        %segment it into wished trials
        cfg = [];
        cfg.length = 2;
        cfg.overlap = 0.5;
        source_time_seg = ft_redefinetrial(cfg,source_time);
        
        if ~exist([mpath,'/source/time/',subjects{i}],'dir')
            mkdir([mpath,'/source/time/',subjects{i}]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%___ Source Power Computation ___%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        if ~exist([mpath,'/source/power/',subjects{i}],'dir')
            mkdir([mpath,'/source/power/',subjects{i}]);
        end
        
        %transfer to frequency domain & compute csd
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'fourier';
        cfg.channel = 'all';
        %cfg.taper = 'hanning';
        cfg.taper = 'dpss';
        cfg.foilim = [1 40];
        cfg.tapsmofrq = 2;
        source_freq_seg = ft_freqanalysis(cfg,source_time_seg);
        
        %get power
        cfg = [];
        psd_for_fooof = ft_freqdescriptives(cfg,source_freq_seg);
        %clean the structure
        psd_for_fooof = rmfield(psd_for_fooof,{'dimord','cumtapcnt','cumsumcnt','cfg'});
        
        %save psd
        save([mpath,'/source/power/',subjects{i},'/',subjects{i},'_source_pow.mat'],'psd_for_fooof');
        
        clear source_time_seg source_time_tmp label
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%___ Parcel Time Computation ___%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Create a new source_time structure with continuous data to compute parcel time course (the previous was packed into subfields calles 'source_xx', which is not suitable here)
        
        %line trials up (only for parcellation -> then again back to trial structure)
        trial_source = horzcat(trial_source{:});
        
        %save results in a structure
        source_time_tmp.pos = parcel.pos;           % we need to go back to the template_grid and its .pos, so we use this here. source_time_tmp.pos WOULD BE the the subject-specific grid (so what is in leadfield)
        source_time_tmp.avg.mom = trial_source;     % #sources X [sampling_rate*trials]
        source_time_tmp.unit = 'cm';
                
        cfg = [];
        cfg.method = 'eig';
        cfg.parcellation = 'mask';              %fieldname with wished parcellation
        cfg.parameter = 'mom';                  %fieldname with data that should be parcellated
        parcel_time = ft_sourceparcellate(cfg,source_time_tmp,parcel);
        
        %check for rank deficiency of parcels
        if rank(parcel_time.mom) ~= size(parcel_time.mom,1); error('Rank deficiency: Parcel time courses are linearly dependent.'); end
        
        %leakage reduction [https://ohba-analysis.github.io/osl-docs/matlab/osl_example_roinets_1_synthetic.html]
        parcel_time.mom = ROInets.remove_source_leakage(parcel_time.mom,'closest');
        parcel_time.mom = mat2cell( parcel_time.mom,size(parcel_time.mom,1),cellfun('size',source_time.trial,2) ); % <- back to original splitting of data (as in source_time or data_select_clean)       
        
        %put it into a structure to segment parcel time courses
        parcel_time.trial = parcel_time.mom;
        parcel_time.time = source_time.time;
        
        %segment it (back) into wished trials (as in line 119, so it matches with source data)
        cfg = [];
        cfg.length = 2;
        cfg.overlap = 0.5;
        parcel_time_seg = ft_redefinetrial(cfg,parcel_time);
        
        %add fields
        parcel_time_seg.brainordinate = parcel_time.brainordinate;
        
        %just change the name
        parcel_time = parcel_time_seg;
        
        %save subject parcel time-courses
        if ~exist([mpath,'/parcel/time/',subjects{i}],'dir')
            mkdir([mpath,'/parcel/time/',subjects{i},'/']);
        end
        
        save([mpath,'/parcel/time/',subjects{i},'/',subjects{i},'_parcel_time.mat'],'parcel_time');
        
    catch
        catch_sub = horzcat(catch_sub,i);
    end
end

warning(['problematic subjects with index: ',num2str(catch_sub)])
