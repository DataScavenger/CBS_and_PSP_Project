%% forward model script: Construct a head model, a 3D grid of cortical sources

%clean workspace and command window
clearvars;
clc;

% path settings
mpath = 'C:/data';                                     %mainpath
ft_path = 'C:/toolboxes/fieldtrip-20201214';           %fieltrip path
fct_path = [mpath,'/functions'];                       %function path (my own functions)
scp_path = [mpath,'/scripts'];                         %script path
cbs_info_path = [mpath,'/cbs_info.mat'];               %cbs_info

%define path to fieldtrip & functions & raw data
addpath(ft_path,fct_path,scp_path);
ft_defaults;

%load info file
load(cbs_info_path)
subjects = fieldnames(cbs_info);

%put 'already_exists' to 1 to skip the calculation of the headmodel from the mri, and only calculate the subjective grid on the template grid.
already_exists = 0;

%choose subject
i = find(strcmpi(subjects,'pd65'));

problems = [];

for i = 1:length(subjects)
    
    close all
    
    try
        
        %path to MRI
        if ~isempty(cbs_info.(subjects{i}).mripath)
            mri_path = cbs_info.(subjects{i}).mripath;
        end
        
        %path to raw meg data
        raw_data_path = [cbs_info.(subjects{i}).rest.path,'.fif']; if iscell(raw_data_path); raw_data_path = [raw_data_path{1},'.fif']; end
        %path to fiducials
        fiducials_path = cbs_info.(subjects{i}).fiducials;
        
        %where would you like to save the results
        if ~exist(['D:/more_clean_data/headmodels/',subjects{i}],'dir')
            mkdir(['D:/more_clean_data/headmodels/',subjects{i}]);
        end
        
        %only enter the hdm calculation if it does not already exist
        if ~exist('already_exists','var') ||  already_exists ~= 1
            
            %read mri
            clear mri
            mri  = ft_read_mri(mri_path);
            mri.anatomy = double(mri.anatomy);
            %raw data file header
            hdr = ft_read_header(raw_data_path);
            %read fiducials
            if ~isempty(fiducials_path)
                fid_ID = fopen(fiducials_path);
                fiducials = textscan(fid_ID,'%f','HeaderLines',2);
                fiducials = reshape(fiducials{1},[3,3])';          %Position in rows: RPA , Nasion, LPA
                fclose(fid_ID);
            end
            
            if strcmpi(subjects{i},'cbs06') || strcmpi(subjects{i},'cbs14')
                %first properly center the MRI (otherwise for these specific subjects the anatomical data are too far off from the origin after ft_volumerealign)
                cfg = [];
                cfg.method = 'interactive';
                cfg.coordsys = 'acpc';      %find anterior commisure and posterior commisure and a positive point on z sagital plane
                mri = ft_volumerealign(cfg,mri);
                %account for the dimension missmatch -> make a [256 256 256] cube
                cfg = [];
                mri = ft_volumereslice([],mri);
            end
            %use volume_reslice for pd67 / or also if you can't easily locate anatomical landmarks
            if strcmpi(subjects{i},'pd67')
                cfg = [];
                mri = ft_volumereslice([],mri);
            end
            
            %align MRI to neuromag coordinate system
            cfg = [];
            cfg.coordsys = 'neuromag';
            if ~isempty(fiducials_path)
                cfg.fiducial.rpa = fiducials(1,:);
                cfg.fiducial.nas = fiducials(2,:);
                cfg.fiducial.lpa = fiducials(3,:);
                cfg.method = 'fiducial';
            else
                cfg.method = 'interactive';
            end
            real_mri = ft_volumerealign(cfg,mri);   %the function does not change the anatomical data but calculates a tranformation matrix to get the anatomical data into the 'MEG'-coordinate system given by the landmarks
            real_mri = ft_convert_units(real_mri,'cm');
            
            %to view the transformed anatomical T1 scan: ft_sourceplot([],real_mri)
            
            %segment the brain into gray, white and csf matter
            cfg = [];
            cfg.write = 'no';
            cfg.keepintermediate = 'no';
            segmented_mri = ft_volumesegment(cfg, real_mri);
            
            %check that segmentation was successful. Important, as from this the headmodel is calculated. We use the (to the MEG-landmarks-coordinate system) aligned MRI, which is important to have a proper placement of the headmodel/brain potate under the MEG sensors
            test = segmented_mri;
            test.anatomy = mri.anatomy;
            cfg = [];
            cfg.funparameter = 'gray';
            ft_sourceplot(cfg,test);
            clear test
            
            %create single shell, realistic headmodel
            cfg = [];
            cfg.method = 'singleshell';
            hdm = ft_prepare_headmodel(cfg,segmented_mri);
        else
            
            %load headmodel & grid (subject specific) to recalculate grid
            load(['D:/more_clean_data/headmodels/',subjects{i},'/',subjects{i},'_forward_model.mat']);
            clear grid leadfield
            
        end
        
        %use Jans cortical surface grid
        load 'C:/data/cortical_grid.mat'
        template_grid.pos = cortical_grid' .* 100; %units are initially in m, now in cm
        template_grid.unit = 'cm';
        clear cortical grid
        
        %create subject grid
        cfg                         = [];
        cfg.method                  = 'basedonmni';  %for wrapping of the grid around the real_mri
        cfg.template                = template_grid; %the template grid
        cfg.sourcemodel.nonlinear   = 'yes';
        cfg.mri                     = real_mri;
        cfg.sourcemodel.unit        ='cm';
        grid                        = ft_prepare_sourcemodel(cfg); %creates a individualized grid per subject
        
        %make a figure of the single subject headmodel, and grid positions
        sens = ft_read_sens(raw_data_path,'senstype','meg');
        if strcmpi(subjects{i},'cbs01'); load('C:/data/new_grad_cbs01.mat'); sens = cbs01_grad_new; end
        if strcmpi(subjects{i},'pd59'); load('C:/data/new_grad_pd59.mat'); sens = pd59_grad_new; end
        
        figure;
        hold on;
        ft_plot_vol(hdm, 'edgecolor', 'none','facecolor','skin', 'facealpha', 0.5);
        ft_plot_mesh(grid.pos(grid.inside,:));
        plot3(grid.pos(:,1),grid.pos(:,2),grid.pos(:,3),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8)
        hold off
        
        figure; hold on;
        ft_plot_vol(hdm);
        ft_plot_sens(sens, 'style', 'b*');
        saveas(gcf,['D:/more_clean_data/headmodels/',subjects{i},'/',subjects{i},'_sensor_hdm_alignment'])
        hold off
        
        %save variables
        cd(['D:/more_clean_data/headmodels/',subjects{i},'/'])
        save([subjects{i},'_forward_model.mat'],'hdm','grid','real_mri');
        
    catch
        problems = horzcat(problems,i);
    end
    
end



