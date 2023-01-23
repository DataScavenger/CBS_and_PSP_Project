%% AAL Parcellation of template grid according to atlas

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

%load cortical grid (from Jan)
load([mpath,'\cortical_grid.mat'])
template_grid.pos = cortical_grid' .* 100; %units: m -> cm
template_grid.inside = ones(length(template_grid.pos),1) == 1;
template_grid.mask = ones(length(template_grid.pos),1) == 1;
template_grid.coordsys = 'mni';
template_grid.unit = 'cm';
clear cortical grid

%load atlas
atlas = ft_read_atlas([ft_path,'/template/atlas/aal/ROI_MNI_V4.nii']);
atlas = ft_convert_units(atlas,'cm');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ___Parcellation___ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove parcels ( -> insula, cingulate gyrus, fusiform gyrus, parahippocampas gyrus, amygdala, ... should be deleted, as the structures are further away from the superficial cortex and no if any grid points are refering to it )
delete_parcels = {'Rectus','Olfactory','Insula','Cingulum','Hippocampus','ParaHippocampal','Amygdala','Thalamus'...
                  'Pallidum','Putamen','Caudate','Vermis','Fusiform','Frontal_Med_Orb_L','Frontal_Med_Orb_R'};
area_labels = atlas.tissuelabel( ~(contains(atlas.tissuelabel,delete_parcels) ));

%sum parcels having the 'characteristics' mentioned in a cell, together
keywords = {{'Precentral_L','Paracentral_Lobule_L'},...     %Precentral_L (note: A cortical sheet is used. The Paracentral Lobule lays medial. Only very few grid points are attributed to it. Thus, they are attributed to the next closet area)
            {'Precentral_R','Paracentral_Lobule_R'},...     %Precentral_R (note: logic as before)
            {'Frontal_Sup_L','Frontal_Sup_Medial_L'},...    %Frontal_Sup_L
            {'Frontal_Sup_R','Frontal_Sup_Medial_R'},...    %Frontal_Sup_R
            {'Parietal_Sup_L','Precuneus_L'},...            %Parietal_Sup_L (note: logic as before)
            {'Parietal_Sup_R','Precuneus_R'},...            %Parietal_Sup_R (note: logic as before)
            {'Heschl_L','Temporal_Sup_L'},...               %Temporal_Sup_L
            {'Heschl_R','Temporal_Sup_R'},...               %Temporal_Sup_R
            {'Occipital_Sup_L','Cuneus_L','Calcarine_L'},...%Occipital_Sup_L
            {'Occipital_Sup_R','Cuneus_R','Calcarine_R'},...%Occipital_Sup_R
            {'Occipital_Inf_L','Lingual_L'},...             %Occipital_Inf_L
            {'Occipital_Inf_R','Lingual_R'},...             %Occipital_Inf_R
            {'Cerebellum_Crus1_L','Cerebellum_Crus2_L','Cerebellum_3_L','Cerebellum_4_5_L','Cerebellum_6_L','Cerebellum_7b_L','Cerebellum_8_L','Cerebellum_9_L','Cerebellum_10_L'},...  %Cerebellum_L
            {'Cerebellum_Crus1_R','Cerebellum_Crus2_R','Cerebellum_3_R','Cerebellum_4_5_R','Cerebellum_6_R','Cerebellum_7b_R','Cerebellum_8_R','Cerebellum_9_R','Cerebellum_10_R'}};    %Cerebellum_R
groupname = {
            'Frontal_Sup_Orb_L','Frontal_Sup_Orb_R',...
            'Frontal_Mid_L','Frontal_Mid_R',...
            'Frontal_Mid_Orb_L','Frontal_Mid_Orb_R',...
            'Frontal_Inf_Oper_L','Frontal_Inf_Oper_R',...
            'Frontal_Inf_Tri_L','Frontal_Inf_Tri_R',...
            'Frontal_Inf_Orb_L','Frontal_Inf_Orb_R',...
            'Rolandic_Oper_L','Rolandic_Oper_R',...
            'Supp_Motor_Area_L','Supp_Motor_Area_R',...
            'Postcentral_L','Postcentral_R',...
            'Parietal_Inf_L','Parietal_Inf_R',...
            'SupraMarginal_L','SupraMarginal_R',...
            'Angular_L','Angular_R',...
            'Temporal_Pole_Sup_L','Temporal_Pole_Sup_R',...
            'Temporal_Mid_L','Temporal_Mid_R','Temporal_Pole_Mid_L','Temporal_Pole_Mid_R',...
            'Temporal_Inf_L','Temporal_Inf_R',...
            'Occipital_Mid_L','Occipital_Mid_R',...
            'Precentral_L',...
            'Precentral_R',...
            'Frontal_Sup_L',...
            'Frontal_Sup_R',...
            'Parietal_Sup_L',...
            'Parietal_Sup_R',...
            'Temporal_Sup_L',...
            'Temporal_Sup_R',...
            'Occipital_Sup_L',...
            'Occipital_Sup_R',...
            'Occipital_Inf_L',...
            'Occipital_Inf_R',...
            'Cerebellum_L',...
            'Cerebellum_R'};

%The updated new area 'groups' (-> #result = #[area_labels - delete_parcels] - #[elements keywords] + #[groups keywords])
area_labels = cbs_sumareas(area_labels,keywords);

%'group_labels' is later used to connect the mask indices with an area 'name'
group_labels = area_labels;
nested_cells = cellfun(@iscell,area_labels);
group_labels(nested_cells) = groupname(nested_cells); %These should be the names behind 'keywords'

%cfg in template grid
cfg = [];
cfg.atlas = atlas;
cfg.maskparameter = 'mask';
cfg.minqueryrange = 1;
cfg.maxqueryrange = 7;

idx_no_areas = [];
idx_multi_areas = [];
pointlabel = cell(length(template_grid.pos),1);

%connect template source with areas
for k = 1:length(template_grid.pos)
    
    %go through points
    point.inside = template_grid.inside(k);
    point.pos = template_grid.pos(k,:); %The position units must match with the point.unit field
    point.mask = template_grid.mask(k);
    point.coordsys = template_grid.coordsys;
    point.unit = 'cm';
    
    %save labels
    mask = ft_volumelookup(cfg,point);
    
    if ~isempty( mask.name(mask.count == 1) )
        
        %save the result
        pointlabel{k} = mask.name(mask.count == 1);
        
        %save multiple assignments
        if length( mask.name(mask.count == 1) ) > 1
            idx_multi_areas = [idx_multi_areas,k];
        end
        %save no label found
        if strcmpi( mask.name(mask.count == 1), 'no_label_found')
            idx_no_areas = [idx_no_areas,k];
        end
    end
    
end
clear mask

%from multiple areas simple choose the first assignment
for k = 1:length(idx_multi_areas)
    pointlabel{idx_multi_areas(k)} = pointlabel{idx_multi_areas(k)}(1);
end
%logical indexing empty cells ('out of range' fieldtrip note)
pointlabel( cellfun(@isempty,pointlabel) ) = {'-'};
%peal off one 'layer' of cell
pointlabel = [pointlabel{:}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ___ PointLabel -> AreaLabel ___%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%connect pointlabel with indices of area labels (depending on in which cell the pointlabel is located)

%make a mask
mask = zeros(length(pointlabel),1);
for k = 1:length(area_labels)
       mask( contains(pointlabel,area_labels{k}) ) = k;
end

%correct grid points (non-mask-attributed points, presumably incorrectly attributed points) -> look up with cbs_parcel_visualization
CorrGrid = {{'Frontal_Sup_L',[554 559]},{'Frontal_Mid_Orb_L',552},...
            {'Frontal_Sup_Orb_L',472},{'Frontal_Sup_Orb_R',[515]},{'Supp_Motor_Area_R',418},...
            {'Frontal_Sup_Orb_R',515},{'Postcentral_L',[251,252]},{'Postcentral_R',[206,248,250]},...
            {'Parietal_Sup_L',[59,176]},{'Occipital_Sup_L',62},{'Occipital_Mid_L',67},...
            {'Occipital_Inf_L',[68,72,536]},{'Occipital_Inf_R',[541,545]}};

for c = 1:length(CorrGrid)
    %extract wished masklabel and corresponding grid points
    region = CorrGrid{c}(1);
    gridpnts = CorrGrid{c}{2};
    %put them into mask
    mask(gridpnts) = find(contains(group_labels,region));
end

%save results in a structure
parcel = [];
parcel.pos = template_grid.pos;
parcel.unit = 'cm';
parcel.coordsys = 'mni';
parcel.mask = mask;                %mask with mask label indices
parcel.masklabel = group_labels';  %area labels

%save the results
save([mpath,'/parcellation.mat'],'parcel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ___Plausibility Plot___ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot to check for plausibility [some region-indications from brainnetome atlas.tissuelabels]

%all areas
check_areas = group_labels;% Or select group labels e.g.: {'Occipital_Med_R','Occipital_Lat_R','Rolandic_Area_L'};
%or a choosen number of areas
check_areas = {'Frontal_Inf_Oper_L','Frontal_Inf_Oper_R','Frontal_Inf_Tri_L','Frontal_Inf_Tri_R'};

%prepare colors
Col = colormap('hsv'); close;
Col = Col(round(linspace(1,size(Col,1),numel(check_areas))), :);

%Load a surface to make a nicer picture
load([ft_path,'/template/anatomy/surface_white_both.mat'])
mesh = ft_convert_units(mesh,'cm');
ft_plot_mesh(mesh,'facealpha',0.1);
view([90,0]);
hold on
for j = 1:numel(check_areas)
    %which areas correspond to the label
    relevant_areas = find(strncmpi(parcel.masklabel,check_areas{j},numel(check_areas{j})));
    %make a mask for a grand region specified in check_areas [sum of related area masks]
    paint = ismember(parcel.mask,relevant_areas);
    
    %plot area
    region = parcel.pos(paint == 1,:); % 10 .*  to transfer from cm to mm ... but we already used ft_convert_units()
    plot3(region(:,1),region(:,2),region(:,3),'o','MarkerFaceColor',Col(j,:),'MarkerEdgeColor',Col(j,:),'MarkerSize',4)
    hold on
end
