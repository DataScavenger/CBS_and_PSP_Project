%% Translate Elements in a vector defined on a Cloud Grid to the opposite Hemisphere
% Input: - labels with the parcel-labels corresponding labels should appear that only differ by '..._L','..._R'
%        - a vector of size n with numbers in that resemble brain activity
% Output: - a vector of brain activity where the activity of one parcel is translated to the correspondingly parcel on the opposite side

function new_vec = cbs_mirror_activity_parcel(labels,vec)
    
    warning('Data must be organized according: Labels X Activity (Rows X Columns)')
    
    %where to put the activity later
    new_vec = vec;
    
    %check left hemispheric labels
    LeLa = labels(endsWith(labels,'_L'));
    
    %exchange activity in 'vec' with the corresponding parcel '..._R'
    for k = 1:length(LeLa)
        %corresponding parcel
        friend = [LeLa{k}(1:end-1),'R'];
        
        %indices in 'vec' for left and right parcel
        iLeft = strcmpi(LeLa{k},labels);
        iRight= strcmpi(friend,labels);
        
        %exchange activity
        ActiLeft = new_vec(iLeft,:);
        ActiRight= new_vec(iRight,:);
        
        new_vec(iRight,:) = ActiLeft;
        new_vec(iLeft,:)  = ActiRight;
        
    end
    
end