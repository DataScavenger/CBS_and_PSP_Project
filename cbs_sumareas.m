%% function to summarize tissue labels to broader regions

%Input 1: area labels / tissue labels. A cell-array. Per cell, containing one string indicating the area
%Input 2: area keyword. A cell array. Summarize to one area per cell, containing one or multiple strings that a area-string of input 1 must contain to be summarized
%Output : new area labels. A cell array, where cells can contain more than one string. Each cell with its related areas indicated by the strings is regarded as one new area

function tmp = cbs_sumareas(area_labels,keywords)
        
        %temporary variable
        tmp = area_labels;
        %save the keyword areas in cell
        items = cell(1,length(keywords));
        %save the keywords without data
        nothing = [];
        %index variable for proper data placement in items
        n = 0;
        %index variable for putting missing data keywords
        m = 0;
        
        for k = 1:length(keywords)
                        
            if any( contains( tmp,keywords{k} ) )
               %update index n
               n = n + 1;
               %save areas related to the keyword in a single cell
               items{n} = tmp( contains( tmp,keywords{k} ) );
               %delete the related areas in the temporary variable
               tmp( contains( tmp,keywords{k} ) ) = [];
            else
               %update index m
               m = m + 1;
               nothing{m} = keywords{k};
            end
            
        end
        
        %save results
        tmp = [tmp,items];
        
        
        if ~isempty(nothing)
           
           nothing_display = nothing{1};
           
           for c = 2:length(nothing)
               nothing_display = [nothing_display,' ',nothing{c}];
           end
                   
           warning(['No atlas areas related to: ',nothing_display])
        
        end
end