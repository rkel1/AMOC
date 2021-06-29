function function_make_AMOC_MAT_file(data_path,name,scen,lat,max_flag)

% function to extract the AMOC timeseries from the given model name,
% scenario, and at the given latitude and save into a .mat file of the
% relevant format. The max_flag attribute dictates whether the maximum of
% the overturning at each depth should be obtained, or if all depths are
% required.
%
% Richard Kelson
% July 2020

% file parameters
home = pwd ;
if contains(home,'aos')
    storage = '/storage3/rkelson/' ;
else
    storage = home ;
end

mfolder = strcat(data_path,name,'/') ;

% script parameters
warning off backtrace
min_years = 498 ; % min length for piControl scenarios (years)
ssp_cutoff = 2100 ; % don't keep any output past this year in SSP scen.

% Check if requested files exist
var_ref_str = 'msft*z' ; % could be msftmz or msftyz
srch_str = sprintf('%s*%s*%s*.nc',var_ref_str,name,scen) ; 
file_list = dir(fullfile(mfolder,srch_str)) ;

% test if .mat file already exists - if so, don't overwrite
if ~isempty(lat)
    sp = sprintf('%s/AMOC %.0f Data/',home,floor(lat)) ;
    save_name = sprintf('%sAMOC_%.0f_%s_%s.mat',sp,floor(lat),scen,name) ;
elseif isempty(lat) && exist('max_flag','var')
    sp = sprintf('%s/AMOC All Data/',home) ;
    if ~max_flag
        save_name = sprintf('%sAMOC_all_field_%s_%s.mat',sp,scen,name) ;
    else
        save_name = sprintf('%sAMOC_all_lats_%s_%s.mat',sp,scen,name) ;
    end
end

if exist(save_name,'file')
    fprintf(1,' Model [%s] Scenario [%s] already exists\n',...
        name,scen)
    return
end

% test if there are any netCDF files from which to extract
if isempty(file_list)
%     warning('No files available for %s, scenario %s',name,scen)
    return
end


% check duration if piControl - needs to be longer than min_years
if strcmpi(scen,'piControl')
    duration = check_duration(file_list) ;
    
    if duration < min_years
        fprintf('Model %s has %.0f years in piControl - requires %d. Skipping\n',...
            name, duration, min_years) ;
        return
    end
end



fprintf(1,'Working on Model [%s] Scenario [%s]\n',name,scen) ;

% create save folder if it doesn't exist
if isempty(lat) && exist('max_flag','var')
    sp = sprintf('%s/AMOC All Data/',home) ; % save path
elseif nargin > 2
    sp = sprintf('%s/AMOC %.0f Data/',home,floor(lat)) ; % save path
elseif nargin > 1
    sp = sprintf('%s/AMOC All Lats Data/',home) ; % save path
end

if ~exist(sp,'dir')
    mkdir(sp)
end

% extract AMOC based on arguments
if ~isempty(lat)
    [t,dt,amoc] = function_extract_amoc(data_path,name,scen,lat) ;
elseif isempty(lat) && exist('max_flag','var')
    [t,dt,amoc,lats,depths] = function_extract_amoc(data_path,name,scen,[],max_flag) ;
    
    % reduce to 100 years if piControl to avoid large file sizes
    if strcmp(scen,'piControl')
        idcs = find(dt <= dt(1) + years(100)) ;
        
        t   = t(idcs) ;
        dt   = dt(idcs) ;
        amoc = amoc(:,:,idcs) ;
    
    end
else
    [t,dt,amoc,lats] = function_extract_amoc(data_path,name,scen) ;
end

% % if piControl and less than min_years years long, discard
% if strcmp(scen,'piControl') && years(range(dt)) < min_years
%     fprintf('Model %s has %.0f years in piControl - requires %d. Skipping\n',...
%         name, years(range(dt)), min_years) ;
%     return
% end

% if SSP 5-85 or 1-26 and past year ssp_cutoff, discard
if strcmp(scen,'ssp585') && any(year(dt) > ssp_cutoff)
    idcs = find(year(dt) < 2100) ;
    
    t = t(idcs) ;
    dt = dt(idcs) ;
    amoc = amoc(idcs) ;
end

% save to the relevant mat file
if ~isempty(lat)
    save_name = sprintf('%sAMOC_%.0f_%s_%s.mat',sp,floor(lat),scen,name) ;
    save(save_name,'t','dt','amoc','name') ; 
elseif isempty(lat) && exist('max_flag','var')
    if ~max_flag
        save_name = sprintf('%sAMOC_all_field_%s_%s.mat',sp,scen,name) ;
        save(save_name,'t','dt','amoc','lats','depths','name') ; 
    else
        save_name = sprintf('%sAMOC_all_lats_%s_%s.mat',sp,scen,name) ;
        save(save_name,'t','dt','amoc','lats','name') ;  
    end
end

end


