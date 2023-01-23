


function A = cbs_cutstring(inp,sym)

A = {};

if ischar(inp)
    
    % where to split the string into segments
    to = strfind(inp,sym);
    from = [1,to + 1];
    % cut segments as defined per column
    index = vertcat(from,[to - 1,length(inp)]);
    
    % index for proper placing in cell array A.
    cou = 1;
    
    for j = 1:size(index,2)
        if ~isempty(inp(index(1,j):index(2,j)))
            
            % Cut out the segment of interest and place it in the cell array
            A{cou} = inp(index(1,j):index(2,j));
            
            cou = cou + 1;
            
        end
    end
    
end


end