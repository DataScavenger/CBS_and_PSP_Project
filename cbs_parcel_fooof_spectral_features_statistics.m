%% Statistical Test on Peaks

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

%load cbs project infos
load([mpath,'/cbs_info.mat']);  %cbs patients info

%load parcellation
load([mpath,'/parcellation.mat']);

%plot statistics
load([ft_path,'/template/anatomy/surface_white_both.mat'])
mesh = ft_convert_units(mesh,'cm');

%load Spectral Features (Peak Oscillation & Amplitude) -> get subject oscillations
load([mpath,'/parcel/ana_shift/spectral_features.mat'],'SpecFeat');

%load data 'CoM'
load([mpath,'/parcel/ana_CoM/CoM.mat']);CoE = tmp; clear tmp %CoM Parcel

%load significant cluster parcels
load([mpath,'/parcel/ana_stat/Significant_Regions.mat']);

%subjects
subjects = SpecFeat.subjects;

%store figures
if ~exist([mpath,'/parcel/ana_stat/'])
    mkdir([mpath,'/parcel/ana_stat/'])
end

%%%%%%%%%%%%%%
%%% Script %%%
%%%%%%%%%%%%%%

%___subject handling

%tremor subjects (from pd)
tremor = cbs_clean_subjects(subjects,cbs_info,{'rest','tremor'},'yes');
%notremor subjects
notremor = cbs_clean_subjects(subjects,cbs_info,{'rest','tremor'},'-yes');

%groups
groups = {'hc',{'cbs','psp'},'pd'};

%group indices
indices.hc = contains(subjects,'hc');
indices.atyp = contains(subjects,{'cbs','psp'});
indices.pd = logical( contains(subjects,'pd') .* contains(subjects,notremor) );

%___parcel handling

%use overlapping clusters of hc vs. atyp & pd vs. atyp
parcels_of_interest = parcel.masklabel( logical( SigR.hcVSaps .* SigR.pdVSaps ) );

%___Oscillations & Amplitude handling

%use amplitude to select major oscillation
use_amp_choices = {'yes','yes'};      %use: {'no','no'} for all peaks + {[4,13],[13,30]}, {'yes','yes'} for the largest peaks + {[4,13],[13,30]}, and {'no'} + {[4,30]}
%choose frequency bin interval for oscillations
intervall_choices = {[4,13],[13,30]};
%store results
StoreLMM = [];

for ch = 1:length(intervall_choices)
    
    intervall = intervall_choices{ch};
    use_amp = use_amp_choices{ch};
    
    %individual features (each row is another subject and there are all the subject related oscillation in)
    IndOsc = cell(length(subjects),length(parcels_of_interest));
    IndAmp = cell(length(subjects),length(parcels_of_interest));
    
    %prune data to parcels of interest & interval
    idx = ismember(SpecFeat.labels,parcels_of_interest);
    for a = 1:length(IndOsc)
        %Spectral Features
        osc_candidates = SpecFeat.Osc{a}(idx);
        amp_candidates = SpecFeat.Amp{a}(idx);
        
        %Check candidates for interval
        for k = 1:length(osc_candidates)
            %Oscillations that fit the intervall criterion
            ChooseThese = logical(( osc_candidates{k} >= intervall(1) ) .* ( osc_candidates{k} <= intervall(2) ));
            %Choose the relevant Oscillations
            if ~isempty(ChooseThese)
                
                switch use_amp
                    
                    case 'yes'
                        
                        [amp,take] = max( amp_candidates{k} .* ChooseThese );  %largest amplitude
                        
                        if amp ~= 0
                            IndOsc{a,k} = osc_candidates{k}(take);             %related center frequency / oscillation
                            IndAmp{a,k} = amp;                                 %related amplitude
                        end
                        
                        
                        
                    otherwise
                        
                        IndOsc{a,k} = osc_candidates{k}(ChooseThese);
                        IndAmp{a,k} = amp_candidates{k}(ChooseThese);
                end
                
            else
                IndOsc{a,k} = [];
                IndAmp{a,k} = [];
            end
            
        end
        
    end
            
    %___Proportion of High / Low frequencies___%
    if intervall(1)==4 && intervall(2)==30
        
        %Beta-Borders (where begins high, where begins los)
        border = 20;
        
        high2low = zeros(size(IndOsc,1),length(border));
        %I simply calculate the ratio (Osc > border(n)) / (all Osc) in the whole ROI
        for n = 1:length(border)
            for k = 1:size(IndOsc,1)
                current_sub_Osc = vertcat( IndOsc{k,:} );
                high2low(k,n) = numel( current_sub_Osc( current_sub_Osc > border(n) ) ) / numel( current_sub_Osc );
            end
        end
    end
        
    %___CoE handling
    
    %parcels of interest indices in the CoE-structure
    idx = ismember(CoE.labels,parcels_of_interest);
    
    %security check if orientation of labels is the same in the CoE as for the parcels_of_interest
    if ~all( strcmpi( CoE.labels(idx),parcels_of_interest ) )
        error('The orientation of labels in parcles_of_interest and CoE.labels is probably different. They must match, so taht the correct CoE values are attributed to the corresponding parcels.');
    end
    
    %get it from the CoE
    CoE_Region = CoE.CoM(idx,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Linear Mixed Model %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %___data table
    
    A = [];
    age = [];
    test = [];
    for s = 1:length(subjects)
        
        %get subject s
        current_subject = subjects(s);
        
        %check if subject is no tremor
        if ~ismember(current_subject,notremor)
            
            warning(['Skip ',current_subject{:},': Not part of notremor group.']);
            
            age(s) = nan;
            %test(s) = nan;
            
            continue
        end
        
        %get subject oscillations & amplitudes
        current_osc = IndOsc(s,:);
        current_amp = IndAmp(s,:);
        
        %group name
        if any( ismember(subjects,current_subject) .* indices.atyp )
            Gr = {'Atyp'};
        end
        if any( ismember(subjects,current_subject) .* indices.hc )
            Gr = {'Hc'};
        end
        if any( ismember(subjects,current_subject) .* indices.pd )
            Gr = {'Pd'};
        end
        
        %get ag (for some Hc the Age was directly entered in the structure under '.age')
        if contains( 'age',fieldnames(cbs_info.(current_subject{:})) )
            
            current_age = cbs_info.(current_subject{:}).age;
            
            if isstr(current_age)
                current_age = str2num(current_age);
            end
            
        else
            
            %birthday and recording date
            birthday = cbs_info.(current_subject{:}).birth;
            recording_day = cbs_info.(current_subject{:}).date;
            
            %get updated cbs_info file <- then delete this if-statement
            if isempty(recording_day)
                catch_sub = vertcat(catch_sub,current_subject{:});
                recording_day = 'xx.xx.2021';
            end
            
            if isempty(birthday)
                catch_sub = vertcat(catch_sub,current_subject{:});
                continue
            end
            
            %get age
            current_age = cbs_GetAge(birthday,recording_day);
            
        end
        
        %extend 'age'
        age(s) = current_age;
        
        for p = 1:length(parcels_of_interest)
            %get subject parcel oscillations & amplitudes
            current_parcel_osc = current_osc{p};
            current_parcel_amp = current_amp{p};
            
            %number of oscillations for parcel{p}
            current_parcel_N = numel(current_parcel_osc);
            
            if isempty(current_parcel_osc)
                current_parcel_osc = nan;
                current_parcel_amp = nan;
                current_parcel_N   = nan;
            end
            
            %get related CoE
            current_parcel_CoE = CoE_Region(p,s);
            
            %information to be added
            B = [repmat(current_subject,length(current_parcel_osc),1),...
                repmat(Gr,length(current_parcel_osc),1),...
                repmat(parcels_of_interest(p),length(current_parcel_osc),1),...
                repmat({current_age},length(current_parcel_osc),1),...
                repmat({current_parcel_CoE},length(current_parcel_osc),1),...
                num2cell( (1:current_parcel_N)' ),...
                num2cell(current_parcel_osc),...
                num2cell(current_parcel_amp),...
                ];
            
            %update A -> later becomes table
            A = vertcat(A,B);
            
            clear B
        end
        
    end
    
    %A = cell2table(A,'VariableNames',{'subject','group','parcel','test','age','CoE','NOsc','Osc','Amp'});
    A = cell2table(A,'VariableNames',{'subject','group','parcel','age','CoE','NOsc','Osc','Amp'});
    
    %make group categorical
    A.group = categorical(A.group);
    A.parcel = categorical(A.parcel);
    A.subject = categorical(A.subject);
    
    A.age = ( A.age  - mean(A.age,'omitnan') );
    A.NOsc= ( A.NOsc - mean(A.NOsc,'omitnan') );
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Response Variable Osc / Amp: LMM with HC, Atyp. PD, Idio. PD %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %Linear Mixed Models
    LMM_Osc_age = fitlme(A,'Osc ~ 1 + group + age + group*age + (1|subject) + (1|subject:parcel) + (1|subject:parcel:NOsc)');   % fixed effects (group intercepts) + random effects (group intercept and age intercept and parcel intercept)
    LMM_Amp_age = fitlme(A,'Amp ~ 1 + group + age + group*age + (1|subject) + (1|subject:parcel) + (1|subject:parcel:NOsc)');   % fixed effects (group intercepts) + random effects (group intercept and age intercept and parcel intercept)
    
    tmp_res = horzcat(LMM_Osc_age.Coefficients.pValue,LMM_Amp_age.Coefficients.pValue);
    
    StoreLMM = horzcat(StoreLMM,tmp_res);
    
    clear tmp_res IndOsc IndAmp
    
end

if length(use_amp_choices) == 2
    %apply multiple comparison correction
    for r = 1:size(StoreLMM,1)
        
        StoreLMM(r,:) = mafdr(StoreLMM(r,:),'BHFDR',true);
        
    end
    
    tmp = mat2cell(StoreLMM,ones(size(StoreLMM,1),1),ones(size(StoreLMM,2),1));
    tmp_row = horzcat( 'exp. variables', LMM_Osc_age.CoefficientNames )';
    tmp_col = { [replace(num2str(intervall_choices{1}),' ','_'),' ',use_amp_choices{1},'_amp',' ','Osc.'],...
        [replace(num2str(intervall_choices{1}),' ','_'),' ',use_amp_choices{1},'_amp',' ','Amp.'],...
        [replace(num2str(intervall_choices{2}),' ','_'),' ',use_amp_choices{2},'_amp',' ','Osc.'],...
        [replace(num2str(intervall_choices{2}),' ','_'),' ',use_amp_choices{2},'_amp',' ','Amp.']};
    
    disp( horzcat( tmp_row, vertcat(tmp_col,tmp) ) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low Freq. Osc / High Freq. Osc Proportion %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if intervall(1)==4 && intervall(2)==30
    
    col = cbrewer('qual', 'Set3', 12, 'pchip');
    
    %center age
    age = age - mean(age,'omitnan');
    
    p_hc = []; p_pd = [];
    for k = 1:size(high2low,2)
        
        high2low_hc = high2low(indices.hc,k);
        high2low_pd = high2low(indices.pd,k);
        high2low_atyp = high2low(indices.atyp,k);
        
        groups = vertcat(repmat({'atyp'},length(high2low_atyp),1),repmat({'hc'},length(high2low_hc),1),repmat({'pd'},length(high2low_pd),1));
        groups = categorical(groups);
        
        age_ = vertcat(age(indices.atyp)',age(indices.hc)',age(indices.pd)');
        
        %make a table for linear model
        C = table( vertcat(high2low_atyp,high2low_hc,high2low_pd),groups,age_,'VariableNames',{'hBetaProp','group','age'} );

        %no need to control for subject or parcel via random effects
        mdl = fitlm(C,'hBetaProp ~ 1 + group + age + group*age');
        
        %save p-value
        p_hc(k) = mdl.Coefficients.pValue(2);
        p_pd(k) = mdl.Coefficients.pValue(3);
    end
    
    %High Beta versus All Oscillations
    disp(['P-Value Proportion High Beta (> 20Hz) compared to all Oscillations: Atyp vs. PD: ',num2str(p_pd),' | Atyp vs. HC: ',num2str(p_hc)])
    
    %Bonferroni Correction
    p_bon = p_pd * length(p_pd);
    p_bon_hc = p_hc * length(p_hc);
    
    pd_border = border( min(find(p_bon < 0.05)) );
    hc_border = border( min(find(p_bon_hc < 0.05)) );
    
    %Where is group distinction prevalent? -> Only enters code if survives bonferroni correction
    if any(p_bon_hc < 0.05)
        
        eventually_these = p_bon_hc .* (p_bon_hc < 0.05);
        eventually_these(eventually_these == 0) = 100;
        
        [~,idx_min] = min(eventually_these);
        
        beta_border = border(idx_min);
        
        jit = -0.175:0.001:0.175;
        
        high2low_hc = high2low(indices.hc,idx_min);
        high2low_pd = high2low(indices.pd,idx_min);
        high2low_atyp = high2low(indices.atyp,idx_min);
        
        figure('Position',[25 25 700 600])
        h = boxplot(vertcat(high2low_hc,high2low_pd,high2low_atyp),...
            vertcat(repmat(1,sum(indices.hc),1),2*repmat(1,sum(indices.pd),1),3*repmat(1,sum(indices.atyp),1)),...
            'Color',[col(1,:);col(6,:);col(4,:)],...
            'symbol','');
        set(h,'linew',1)
        hold on;
        scatter(ones(length(high2low_hc),1) + jit(1,randperm(length(jit),sum(indices.hc)))',high2low_hc,50,'MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:)); hold on;
        scatter(2*ones(length(high2low_pd),1) + jit(1,randperm(length(jit),sum(indices.pd)))',high2low_pd,50,'MarkerFaceColor',col(6,:),'MarkerEdgeColor',col(6,:)); hold on;
        scatter(3*ones(length(high2low_atyp),1) + jit(1,randperm(length(jit),sum(indices.atyp)))',high2low_atyp,50,'MarkerFaceColor',col(4,:),'MarkerEdgeColor',col(4,:)); hold on;
        title('')
        ylabel('Beta Proportion [>20Hz]')
        ylim([-0.09,1])
        drawbrace([1 0.925],[3 0.925],5,'Color',[0 0 0],'LineWidth',2); text([1.933 1.933],[0.95 0.95],'**','FontSize',20);
        drawbrace([2 0.825],[3 0.825],5,'Color',[0 0 0],'LineWidth',2); text([2.4675 2.4675],[0.85 0.85],'*','FontSize',20);
        legend({'HC','PD','APS'},'FontWeight','bold','Location','northeastoutside')
        set(gca,'XTickLabel','')
        set(gca,'FontSize',20)
        
        print(gcf,[mpath,'/parcel/ana_shift/HighBeta_Proportion.tiff'],'-dtiff','-r300')
        
    end
    close all
    
end
