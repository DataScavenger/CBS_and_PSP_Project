%function to clean the 'subjects' cell-array that I frequently use.

%___input-variables
%   subjects := cell array that contains a subject ID in each cell (character array)
%   cbs_info := info structure with all subjects as fieldnames. Behind each subject field
%               there should be information that helps the function to clean the
%               subjects array accordingly.
%   ___name-value pairs
%      'some fieldname','characteristic'    :   E.g.: ...'eyes','open',... to choose all subjects having the eyes open
%      'some fieldname','-characteristic'   :   E.g.: ...'tremor','-yes',... to exclude all subjects having tremor
%      {'some fieldname','some fieldname'},'characteristic' : reach nested fieldnames and their properties. E.g.: ... {'rest','eyes'},'open' to include these subjects
%      {'some fieldname','some fieldname'},'-characteristic': reach nested fieldnames and their properties. E.g.: ... {'rest','eyes'},'-open' to exclude these subjects
%      'exclude',{'s01','s02'}              :   List of subjects that should be included in any case
%      'include',{'s04','s11'}              :   List of subjects that should be included in any case
%
% Note: Be sure that the fieldnames do exist for all subjects. Otherwise the function will crash.

function [S,tmp] = cbs_clean_subjects(subjects,cbs_info,varargin)

%output variables
tmp = zeros( length(subjects), length(varargin)/2); %tmp: subjects X |name-value-pairs|
S = [];

%check fieldnames of interest -> fields loop ______________________________
for j = 1:length(varargin)/2
    
    %check variable input 'name' j
    fld = varargin{2*j-1};
    %check if it is character array. If so -> make {}
    if ischar(fld)
        fld = {fld};
    end
    
    switch fld{1}
        case 'exclude' %exclude no matter what other characteristics
            tmp( ismember(subjects,varargin{2*j}), j ) = -inf;
            tmp( ~ismember(subjects,varargin{2*j}), j ) = 1;
        case 'include' %include no matter what other characteristics
            tmp( ismember(subjects,varargin{2*j}), j ) = inf;
            tmp( ~ismember(subjects,varargin{2*j}), j ) = 1;
        otherwise
            %check variable input 'value' j
            chz = varargin{2*j};
            %check if it is a character array
            if ~ischar(chz)
                error('name-value-pair: value must be datatype "character" ')
            end
            
            %check characteristic (per subject) -> subject loop ___________
            for i = 1:length(subjects)
                %check field / nested field (per subject)
                a = getfield(cbs_info.(subjects{i}),fld{:});
                %compare with characteristic
                if startsWith(chz,'-')
                    if contains(a,chz(2:end))
                        tmp(i,j) = 0;
                    else
                        tmp(i,j) = 1;
                    end
                else
                    if contains(a,chz)
                        tmp(i,j) = 1;
                    else
                        tmp(i,j) = 0;
                    end
                end
            end
    end
end

% summarize all information in tmp
S = subjects( sum(tmp,2) >= size(tmp,2) );

end