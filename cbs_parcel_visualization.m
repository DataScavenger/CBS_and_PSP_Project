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

%load source information (for positions, unit, labels)
load([mpath,'/parcellation.mat']);

%plot mesh
load([ft_path,'/template/anatomy/surface_white_both.mat'])
mesh = ft_convert_units(mesh,'cm');

%%%%%%%%%%%%%%%%%%%%%%
%%% Color Regions %%%%
%%%%%%%%%%%%%%%%%%%%%%

%regions
parcels_of_interest = {{'Rolandic_Oper_L','Rolandic_Oper_R',...
                        'Postcentral_L','Postcentral_R','Precentral_L','Precentral_R'},...
                       {'Parietal_Sup_L','Parietal_Sup_R',...
                        'Parietal_Inf_L','Parietal_Inf_R',...
                        'SupraMarginal_L','SupraMarginal_R',...
                        'Angular_L','Angular_R'},...
                       {'Frontal_Sup_Orb_L','Frontal_Sup_Orb_R','Frontal_Mid_L','Frontal_Mid_R',...
                        'Frontal_Mid_Orb_L','Frontal_Mid_Orb_R','Frontal_Inf_Oper_L','Frontal_Inf_Oper_R',...
                        'Frontal_Inf_Tri_L','Frontal_Inf_Tri_R','Frontal_Inf_Orb_L','Frontal_Inf_Orb_R',...
                        'Supp_Motor_Area_L','Supp_Motor_Area_R',...
                        'Frontal_Sup_L','Frontal_Sup_R'},...
                       {'Temporal_Pole_Sup_L','Temporal_Pole_Sup_R','Temporal_Mid_L','Temporal_Mid_R',...
                        'Temporal_Pole_Mid_L','Temporal_Pole_Mid_R','Temporal_Inf_L','Temporal_Inf_R',...
                        'Temporal_Sup_L','Temporal_Sup_R'},...
                       {'Occipital_Mid_L','Occipital_Mid_R','Occipital_Sup_L','Occipital_Sup_R',...
                        'Occipital_Inf_L','Occipital_Inf_R'}};

parcels_of_interest = parcel.masklabel;

%tmp
tmp.pos = parcel.pos;
tmp.unit = parcel.unit;
tmp.region = zeros(length(parcel.mask),1);

for i = 1:length(parcels_of_interest)
    if ~ischar(parcels_of_interest{i})
        for j = 1:length(parcels_of_interest{i})
            
            %identify label of parcel in 'parcel'
            label = parcels_of_interest{i}{j};
            %find mask of label
            mask = find( strcmpi(parcel.masklabel,label) );
            
            %identify mask in .pos
            tmp.region( parcel.mask == mask ) = i;
        end
    else
        tmp.region = parcel.mask;
        break
    end
end

remaining = find(tmp.region == 0);
area_attached = (tmp.region ~= 0);
for i = 1:length(remaining)
    %distance
    distance = sqrt( sum((tmp.pos - tmp.pos(remaining(i),:)) .* (tmp.pos - tmp.pos(remaining(i),:)),2) );
    %finding closest non-zero
    overlap = area_attached .* distance;
    overlap(overlap == 0) = 100;
    %...
    [~,idx] = min( overlap );
    
    %replace 0 with the closest 'area'mask
    tmp.region( remaining(i) ) = tmp.region(idx);
end

%plot regions
Cregion = cbrewer('div','Spectral',64);

%interpolate
cfg = [];
cfg.parameter = 'region';
cfg.downsample = 3;
cfg.method = 'cubic';
tmp_interpol = ft_sourceinterpolate(cfg,tmp,mesh);

%plot regions
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'region';
cfg.projmethod = 'nearest';
cfg.funcolormap = Cregion;
cfg.camlight = 'no';
ft_sourceplot(cfg,tmp_interpol)
colorbar off
set(gcf,'Position',[10 10 800 800])
%change views and make images (right -> left -> top)
view([-90 0])
print([mpath,'/parcel/Parcellation_left.tiff'],'-dtiff','-r300');
%other side ...
view([0 90])
print([mpath,'/parcel/Parcellation_top.tiff'],'-dtiff','-r300');
%save image
view([90 0])
print([mpath,'/parcel/Parcellation_right.tiff'],'-dtiff','-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Color Sources of Parcels %%%% -> Look out for mismatches and eventually manually correct for it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parcels
parcels_of_interest = parcel.masklabel;

%color
col = cbrewer('qual','Set3',length(parcels_of_interest),'pchip');

%no label grid points
nogridpnts = parcel.pos( parcel.mask == 0 ,: );
notxt = find( parcel.mask == 0 );

figure('Position',[50 50 1000 600])
ft_plot_mesh(mesh,'facealpha',0.2)
view([0,90])
hold on
for k = 1:length(parcels_of_interest)
    
    %get parcel label
    label = parcels_of_interest{k};
    
    %get parcel grid points
    gridpnts = parcel.pos( parcel.mask == find(strcmpi(parcel.masklabel,label)) ,:);
    
    %grid indices
    gridtxt = find( parcel.mask == find(strcmpi(parcel.masklabel,label)) );
    
    %put the sensor in the plot
    h1 = scatter3(gridpnts(:,1),gridpnts(:,2),gridpnts(:,3),'MarkerFaceColor',col(k,:),'MarkerEdgeColor',col(k,:));
    hold on
    for l = 1:length(gridtxt)
    ht(l) = text(mean(gridpnts(l,1)),mean(gridpnts(l,2)),mean(gridpnts(l,3)),num2str(gridtxt(l)),'FontSize',8);
    hold on
    end
    
    %put parcel label into the corner
    ht(end + 1) = text(8,8,8,replace(label,'_',' '),'FontSize',12);
    hold on
    
    %add grid points that got attributed to no parcel -> to see if they actually nicely fill in with a parcel
    h2 = scatter3(nogridpnts(:,1),nogridpnts(:,2),nogridpnts(:,3),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0.5 0.5 0.5]);
    hold on
    for l = 1:length(notxt)
    htg(l) = text(mean(nogridpnts(l,1)),mean(nogridpnts(l,2)),mean(nogridpnts(l,3)),[' ',num2str(notxt(l))],'FontSize',8);
    hold on
    end
    
    pause(2.)
    
     delete(ht)
     delete(htg)
     delete(h1)
     delete(h2)
    
end

