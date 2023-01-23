%% function to determine 'source' neighbours on a cortical grid (e.g. Jan's coritcal grid) -> Had not been tested on 3D Grid
%
% Input:

% Cortical Grid aka 'cortical_grid' : -> Need .pos and .labels [position data must match with labels by indices]
% Distance      aka 'grid_distance' : -> Distance (in units of the grid) to determine neighbours
% Third agument forwards 'source'

% OR

% Cortical Grid: Parcel-Structure that contains .pos, .masklabel, .mask -> Will give out all direct neighbours of a parcel
% Distance      aka 'grid_distance' : -> Distance (in units of the grid) to determine neighbours
% Third argument forwards 'parcel'
% (Eventual) fourth argument: 'same' -> Correpsonding areas in both areas have the same corresponding neighbours within their hemisphere

function neigh = cbs_prepare_neighbours(cortical_grid,grid_distance,datatype,varargin)
    
if strcmpi(datatype,'source')
    %neighbourhood structure
    for k = 1:length(cortical_grid.labels)
        
        %___put label name
        neigh(k).label = cortical_grid.labels{k};
        
        %___find neighbours: use 'grid_distance' to determine the neighbourhood
        
        %reference position
        reference = cortical_grid.pos(k,:);
        
        %get neighbours (euclidian distance to reference)
        idx =  sqrt( sum( ( cortical_grid.pos - reference ).^2 ,2) ) <= grid_distance;
        %remove the reference as a potential candidate
        idx(k) = false;
        
        %get neighbours
        neigh(k).neighblabel = {cortical_grid.labels{idx}}';
    end
    
elseif strcmpi(datatype,'parcel')
    
    for k = 1:length(cortical_grid.masklabel)
        
        neigh(k).label = cortical_grid.masklabel{k};
        
        %candidates
        candi = cortical_grid.pos(cortical_grid.mask == k,:);
        %potential neighbours labels
        tmp_neigh_labels = {};
        %potential neighbours position
        tmp_neigh_pos = cortical_grid.pos;
        tmp_neigh_pos(cortical_grid.mask == k,:) = nan;
        
        for p = 1:size(candi,1)
            
            %check distance
            idx = sqrt( sum((tmp_neigh_pos - candi(p,:)).^2,2) ) <= grid_distance;
            
            %only regard grid points that are actually part of a parcel
            neigh_mask = cortical_grid.mask(idx);
            neigh_mask = neigh_mask(neigh_mask ~= 0);
            
            %parcel labels
            tmp_neigh_labels{p} = cortical_grid.masklabel( neigh_mask );
            
        end
        
        %get neighbours
        neigh(k).neighblabel = unique( vertcat(tmp_neigh_labels{:}) );
        
    end

    if ( nargin > 3 ) && strcmpi(varargin{1},'same')
        
        for k = 1:length(neigh)
            
            %get current label
            tmp_label_ref = neigh(k).label;

            %look for match on other hemisphere
            if endsWith(tmp_label_ref,'L')
                %match label
                tmp_label_match = tmp_label_ref; tmp_label_match(end) = 'R';
                %match index
                idx = contains({neigh.label}, tmp_label_match);
                
                %merge their neighbours without the hemisphere letter '_L' or '_R'
                tmp = unique( vertcat( replace(neigh(k).neighblabel,{'_L','_R'},''),replace(neigh(idx).neighblabel,{'_R','_L'},'') ) );
                
                %give each hemisphere what they both have together in sum
                
                %left hemisphere
                tmp_L = horzcat( tmp, repmat({'_L'},numel(tmp),1) );
                neigh(k).neighblabel = join(tmp_L,'');
                
                %right hemisphere
                tmp_R = horzcat( tmp, repmat({'_R'},numel(tmp),1) );
                neigh(idx).neighblabel = join(tmp_R,'');
                
            end
            
        end
        
    else
        
    end
    
else
    error('Choose source or parcel as third argument (-> datatype)')
end

end