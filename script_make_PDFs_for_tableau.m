% Script to generate the variables of different models for Tableau
%
% Richard Kelson
% June 2021

clear all

% set file parameters
home = pwd ;
save_path = sprintf('%s/Tableau/',home) ;
data_path = sprintf('%s/AMOC 26 Data/',home) ;

save_name = strcat(save_path,'PDF.xlsx') ;

% set script parameters
scen = 'piControl' ;

% get model names from model_path directory
models = get_model_names(data_path,scen) ; 

T = 498*12 ;


%% load mean of model data
segments = 12.*[5:5:35] ;  

edges = linspace(-0.5, 0.5, 100) ;

E = length(edges) ;
M = length(models) ;
S = length(segments) ;

% preallocate
dtb_array  = NaN(M,E,S) ;
std_array  = NaN(M,E,S) ;

for m_idx = 1:M
    
    % load this model
    model_name = models{m_idx} ;
    load_path = sprintf('%sAMOC_26_%s_%s.mat',data_path,scen,model_name) ;
    
    load(load_path,'amoc','dt','t') ;
    
    % take first 498 years only
    dt   = dt(1:T) ;
    t    = t(1:T) ;
    amoc = amoc(1:T) ;
    
    % get interannual component for varying annual cycle
    interannual = function_x12_filter(dt,amoc) ;
    
    % get PDF of given segment lengths
    for k_idx = 1:S 
        k = segments(k_idx) ;
        [array,~,dtb] = function_trend_id_pdf(t,dt,interannual,k,edges) ;
        
        % put in array
        dtb_array(m_idx,:,k_idx) = dtb ;
        std_array(m_idx,:,k_idx) = std(array.trends) ;
    end
    
end
    
% get median/mean distribution for each segment length, and mean, std
median_dtb  = squeeze(mean(dtb_array,1,'omitnan')) ;
median_std  = squeeze(mean(std_array,1,'omitnan')) ;

%% save into mat file
save_mat = 'PDF.mat' ;

save(save_mat,'edges','segments','dtb_array','median_std','median_dtb') ;



%% tableau reformatting
home = pwd ;
save_path = sprintf('%s/Tableau/',home) ;
data_path = sprintf('%s/AMOC 26 Data/',home) ;
save_name = strcat(save_path,'PDF.xlsx') ;

models = get_model_names(data_path,'piControl') ; 

load('PDF.mat')

E = length(edges) ;
M = length(models) ;
S = length(segments) ;

% reshape matrix into long-form and label
median_dtb_array = NaN(3,E*S) ;
median_std_array = NaN(3,E*S) ;
dtb_array_long = NaN(4,E*S*M) ;


% add dtb array and median std_array
median_dtb_array(3,:) = reshape(median_dtb,1,E*S) ;
median_std_array(3,:) = reshape(2.*median_std,1,E*S) ; % switch to 2sigma
dtb_array_long(4,:)   = reshape(dtb_array,1, E*S*M) ;

% add dimensions to first two columns
edges_col = reshape(repmat(edges,1,S),1,E*S) ;
seg_col   = reshape(repmat(segments./12,E,1),1,E*S) ; % convert to years

median_dtb_array(1,:) = edges_col ;
median_dtb_array(2,:) = seg_col ;

median_std_array(1,:) = edges_col ;
median_std_array(2,:) = seg_col ;

% and for all-model table
model_col = reshape(repmat(models',1,E*S),1,E*S*M) ;
edges_col = reshape(repmat(edges,M,S),1,E*S*M) ;
seg_col   = reshape(repmat(segments./12,E*M,1),1,E*S*M) ; % convert to years


dtb_array_long(2,:) = edges_col ;
dtb_array_long(3,:) = seg_col ;

dtb_array_long = num2cell(dtb_array_long) ;
dtb_array_long(1,:) = model_col ;

% add header labels
median_PDF = vertcat({'Edges','Segment Length','PDF Value'},num2cell(median_dtb_array')) ;
median_STD = vertcat({'Edges','Segment Length',  'Std Dev'},num2cell(median_std_array')) ;
models_PDF = vertcat({'Models','Edges','Segment Length','PDF Value'},dtb_array_long') ;

% save into table
writecell(median_PDF, save_name,'Sheet',1) ;
writecell(models_PDF, save_name,'Sheet',2) ;
writecell(median_STD, save_name,'Sheet',3) ;




