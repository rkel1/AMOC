function [trend_lengths_dep, model_list, sig_levels_dep, std_of_PDFs] = ...
    function_sig_vs_lgth(VAR_FLAG,R_FLAG,scen,trend_lengths_indep,sig_levels_indep)

% Function to execute script_find_sig_trend for either a varying monthly
% climatology (VAR_FLAG = 1) or a constant monthly climatology, at the
% significance levels provided by interpolating from the input trend
% lengths. R_FLAG sets whether the RAPID trend is extracted using VAC or
% CAC methodology.
%
% Outputs:
%   tnd_lengths: vector of lengths the RAPID trend must persist to be
%                   significant at the corresponding significance levels
%   mean_mdl:     mean over model of tnd_lengths
%
% Richard Kelson
% July 2020

% file parameters
home = pwd ; 
data_path = sprintf('%s/AMOC 26 Data/',home) ;

% get RAPID data
rapid_data = load('RAPID_trend.mat') ;


% script parameters
% scen = 'piControl' ; % comparing to Pre-Industrial control scenarios
edges = 0.6.*linspace(-1,1,100) ; % points at which to calculate the fit
if R_FLAG % load RAPID trend
    rapid_trend = rapid_data.trend_var ;
else
    rapid_trend = rapid_data.trend_cte ;
end

% find if multiple significances requested
SIG_FLAG = 0 ;
if nargin > 3
    SIG_FLAG = 1 ;
end


% get model names from model_path directory
all_models = get_model_names(data_path,scen) ;

%%% debug: don't use E3SM suite %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e3sm_idcs = contains(all_models,'E3SM') ;
models = all_models(~e3sm_idcs) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(scen,'piControl')
    T = 498*12 ; % months in 498 years
elseif strcmpi(scen,'historical')
    T = 160*12 ; % months in 160 years
end

%% Find trend length for each model

% preallocate
sig_levels_dep    = NaN(length(models),length(trend_lengths_indep)) ;
std_of_PDFs       = NaN(length(models),length(trend_lengths_indep)) ;
model_list        = cell(length(models),1) ;
if SIG_FLAG
    trend_lengths_dep = NaN(length(models),length(sig_levels_indep)) ;
end

for m = 1:length(models) % iterate models
    model_name = models{m} ;
    
    fprintf(1,'Now working on model %s\n',model_name) ;
    
    load_path = sprintf('%sAMOC_26_%s_%s.mat',data_path,scen,model_name) ;
    load(load_path,'amoc','dt','t') ;
    
    % add this model to model list
    model_list{m} = model_name ;
    
    % take first 498 years only
    dt   = dt(1:T) ;
    t    = t(1:T) ;
    amoc = amoc(1:T) ;

    
    % preallocate
    sig_vector_dep = NaN.*trend_lengths_indep ; % significance of RAPID
    std_dev_of_PDF = NaN.*trend_lengths_indep ; % std dev of each PDF
    for j = 1:length(trend_lengths_indep) % iterate length of segment
        segment_length = trend_lengths_indep(j) ;
        segment_in_months = segment_length * 12 ;
        
        % find distribution of this segment length
        p_arr = function_trend_id_pdf(t,dt,amoc,segment_in_months,edges) ;
        
        % find std dev of PDF in Sv/yr
        std_dev_of_PDF(j) = std(p_arr.trends) ;
        
        % find significance of RAPID trend in std dev units
        sig = abs(rapid_trend / std_dev_of_PDF(j)) ;
        
        % put into significance vector
        sig_vector_dep(j) = sig ;
    end % j
    
    % save into array
    sig_levels_dep(m,:) = sig_vector_dep ;
    std_of_PDFs(m,:)    = std_dev_of_PDF ;
    
    % get trend length required for multiple significances, if requested
    if SIG_FLAG
        intp_length = interp1(sig_vector_dep,trend_lengths_indep,...
                    sig_levels_indep,'linear','extrap') ;
        trend_lengths_dep(m,:) = intp_length ;
    end
    
end

% keep only models that produces output
good_idcs = any(~isnan(sig_levels_dep),2) ;

sig_levels_dep = sig_levels_dep(good_idcs,:) ;
model_list = model_list(good_idcs) ;
if SIG_FLAG
    trend_lengths_dep = trend_lengths_dep(good_idcs,:) ;
end

end