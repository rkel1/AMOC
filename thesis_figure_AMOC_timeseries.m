% Script to plot the AMOC timeseries. 
%
% a) 2004-2018 observed AMOC decomposition into I,C,R components
% b) 1870-2014 Mean Groups 1-3 + std, Caesar and RAPID I components
%
% Richard Kelson
% January 2021

clear all

fprintf(1,'Running %s:\n%s\n',mfilename,repmat('=',1,70)) ;

%% Setup
% file parameters
home = pwd ;
load_data_path = sprintf('%s/AMOC 26 Data/',home) ;
save_data_path = sprintf('%s/AMOC Grouped Timeseries/',home) ;
save_fig_path = sprintf('%s/Figures/Thesis Main/',home) ;

if ~exist(save_data_path,'dir')
    mkdir(save_data_path)
end
if ~exist(save_fig_path,'dir')
    mkdir(save_fig_path)
end

% see if modelsdata is already made - if not, make it
filename = strcat(save_data_path,'AMOC Grouped Model Timeseries.mat') ;
if isempty(dir(filename))
    make_grouped_timeseries_matrix(load_data_path,save_data_path)
end

% Load grouped models
mdls = load(filename) ;

% Load RAPID Data
r = load('RAPID_trend.mat') ;

% Load Caesar Data
cs = load('Caesar_AMOC.mat') ;

% Load bathymetry data
bt = load(strcat(home,'/Bathymetry Data/bathymetry.mat')) ;

%% Set up Figure
tms_fig = figure('Position',[142.6000 308.2000 947 629]) ;
plt = tiledlayout(2,2,'TileSpacing','None') ;


%% Panel 1: Geographic Location
nexttile(1) ; hold on

% set up mapping axes
ax = axesm('mercator',...
    'MapLatLimit',bt.lats([1,end]),...
    'MapLonLimit',bt.lons([1,end]),...
    'Grid','on') ;

% plot bathymetry
h = geoshow(bt.lats,bt.lons,bt.bt','DisplayType','mesh') ; hold on
shading flat
cmocean('topo','pivot',0)

tightmap on
plabel on
mlabel on
mlabel south

% plot RAPID array
plotm([26.5 26.5],[-80 -71],...
      'r',...
      'LineWidth',2) ;
plotm([26.5 26.5],[-30 -15],...
      'r',...
      'LineWidth',2) ;  
        
% text of RAPID array
textm(27.4,-70,'RAPID 26.5^oN',...
    'VerticalAlignment','middle',...
    'FontSize',12,...
    'Color','r',...
    'FontWeight','bold',...
    'FontName','Arial')

title('Location of the RAPID array')
   
text(.01,.99,'a)',...
    'FontSize',20,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'Units','Normalized',...
    'Color','w')

%% Panel 2: RAPID Timeseries
nexttile(2) ; hold on

% legend labels
xstr = '\makebox[1cm]{\textbf{X(t)}} - Raw' ;
istr = '\makebox[1cm]{\textbf{I(t)}} - Interannual' ;
cstr = '\makebox[1cm]{\textbf{C(t)}} - Annual Cycle' ;
rstr = '\makebox[1cm]{\textbf{R(t)}} - Residual' ;

plot(r.times,r.timeseries,'DisplayName',xstr)
plot(r.r_dt,r.var_inter,'k','LineWidth',2,'DisplayName',istr)
plot(r.r_dt,r.var_annual,'r','LineWidth',2,'DisplayName',cstr) 
plot(r.r_dt,r.var_resid,'Color',[0 .6 0],'DisplayName',rstr) 

legend('location','eastoutside','Interpreter','latex','Box','off')
title('AMOC Strength by frequency component')
xlabel('Date')
ylabel('Strength [Sv]')

% set panel lettering
text(.01,.99,'b)',...
        'FontSize',20,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top',...
        'Units','Normalized',...
        'Color','k')

%% Panel 3: X12 Decomposition
nexttile(3) ; hold on

xistr = '\makebox[1cm]{\textbf{X-12}}' ;
xcstr = '\makebox[1cm]{\textbf{X-12}} - C(t)' ;
cistr = '\makebox[1cm]{\textbf{CAC}}' ;
ccstr = '\makebox[1cm]{\textbf{CAC}} - C(t)' ;

plot(r.r_dt,r.var_inter,'Color',[.2 .2 .6],'DisplayName',xistr)
plot(r.r_dt,r.var_annual,'Color',[.2 .2 .6],'HandleVisibility','off') 

plot(r.r_dt,r.cte_inter,'Color',[.6 0 0],'DisplayName',cistr)
plot(r.r_dt,r.cte_annual,'Color',[.6 0 0],'HandleVisibility','off') 

% add labels for I(t), C(t)
text(datetime(2018,1,1),20,'\textbf{I(t)}',...
     'FontSize',10,...
     'Interpreter','latex',...
     'HorizontalAlignment','right')
text(datetime(2018,1,1),7.5,'\textbf{C(t)}',...
     'FontSize',10,...
     'Interpreter','latex',...
     'HorizontalAlignment','right')

legend('location','west','Interpreter','latex','Box','off')
title('Frequency component differences')
xlabel('Date')
ylabel('Strength [Sv]')

ylim([-7 24])

% set panel lettering
text(.01,.99,'c)',...
        'FontSize',20,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top',...
        'Units','Normalized',...
        'Color','k')
    
    
%% Panel 4: Representative Models
nexttile(4) ; hold on

mdls = get_model_names(strcat(home,'/AMOC 26 Data/'),'piControl') ;
T = 100*12 ;
dt = datetime(0,1,1):calmonths(1):datetime(99,12,31) ;

% label representative models
r_mds = {'ACCESS-CM2','CNRM-ESM2-1','NorESM2-MM'} ;
m_idcs = find(contains(mdls,r_mds)) ;

for i = 1:length(m_idcs)
    srch_str = sprintf('%s/AMOC 26 Data/*piControl_%s.mat',home,r_mds{i}) ;
    file_path = dir(fullfile(srch_str)).name ;
    
    pi_sim = load(sprintf('%s/AMOC 26 Data/%s',home,file_path)) ;
    
    % keep first 100 years only
    dt_sim    = pi_sim.dt(1:T) ;
    t     = pi_sim.t(1:T) ;
    amoc = pi_sim.amoc(1:T) ;
    
    % get interannual component
    inter = function_x12_filter(dt_sim,amoc) ;
    
    % plot
    p = plot(dt,inter,'DisplayName',r_mds{i}) ;
end

xlim([datetime(0,1,1) datetime(100,1,1)])

% Add legend
l = legend('Location','eastoutside','box','off') ;

% Add labels
title('100 years of PI timeseries')
xlabel('Year')
ylabel('AMOC Strength [Sv]')

% set panel lettering
text(.01,.99,'d)',...
        'FontSize',20,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top',...
        'Units','Normalized',...
        'Color','k')
        
%% Save
sp = strcat(save_fig_path,'AMOC Timeseries - Breakdown.png') ;
exportgraphics(tms_fig,sp)
    
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_grouped_timeseries_matrix(load_path,save_path)

    % set script parameters
    scen = 'historical' ;

    % get model names from model_path directory
    all_models = get_model_names(load_path,scen) ;

    %%% debug: don't use E3SM suite %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e3sm_idcs = contains(all_models,'E3SM') ;
    models = all_models(~e3sm_idcs) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load grouped models
    gps = load('Grouped Models.mat') ;

    % set datetime for 1870-2014 panel
    dt = datetime(1870,1,1):calmonths(1):datetime(2014,12,1) ;

    %% Get grouped model timeseries
    G1 = length(gps.grp1_models) ;
    G2 = length(gps.grp2_models) ;
    G3 = length(gps.grp3_models) ;

    T = length(dt) ;
    M = length(models) ;

    % preallocate model group timeseries matrices
    g1 = NaN(G1,T) ;
    g2 = NaN(G2,T) ;
    g3 = NaN(G3,T) ;

    % set counters for insertion into right matrices
    g1_c = 1 ;
    g2_c = 1 ;
    g3_c = 1 ;

    % put models into appropriate group matrix
    for m = 1:M
        model_name = models{m} ;
        fprintf(1,'Now working on model %s\n',model_name) ;

        % load data
        search_str = sprintf('*%s*%s*.mat',scen,model_name) ;
        mat_file = dir(fullfile(load_path,search_str)) ;
        mat_path = sprintf('%s%s',load_path,mat_file(1).name) ;
        mdl = load(mat_path) ;

        % extract interannual trend
        amoc_i = function_x12_filter(mdl.dt,mdl.amoc) ;

        % Trim to required dates (1870 - 2014)
        t_idcs = find(year(mdl.dt) >= year(dt(1))) ;
        amoc_dt = mdl.dt(t_idcs) ;
        amoc_i = amoc_i(t_idcs) ;
        
        L = length(t_idcs) ;
        
        % convert model_name to valid variable name
        var_name = matlab.lang.makeValidName(model_name) ;

        % add to relevant group matrix
        if ismember(model_name,gps.grp1_models)
            g1(g1_c,1:L) = amoc_i ;
            grp1.mdl_timeseries.(var_name).timeseries = amoc_i ;
            grp1.mdl_timeseries.(var_name).times = amoc_dt ;
            
            g1_c = g1_c + 1 ;
            
        elseif ismember(model_name,gps.grp2_models)
            g2(g2_c,1:L) = amoc_i ;
            grp2.mdl_timeseries.(var_name).timeseries = amoc_i ;
            grp2.mdl_timeseries.(var_name).times = amoc_dt ;
            
            g2_c = g2_c + 1 ;
        elseif ismember(model_name,gps.grp3_models)
            g3(g3_c,1:L) = amoc_i ;
            grp3.mdl_timeseries.(var_name).timeseries = amoc_i ;
            grp3.mdl_timeseries.(var_name).times = amoc_dt ;
            
            g3_c = g3_c + 1 ;
        end
    end

    % Calculate the mean across models
    grp1.mean = squeeze(mean(g1,1,'omitnan')) ;
    grp2.mean = squeeze(mean(g2,1,'omitnan')) ;
    grp3.mean = squeeze(mean(g3,1,'omitnan')) ;

    % Calculate the std dev across models
    grp1.std = std(g1,'omitnan') ;
    grp2.std = std(g2,'omitnan') ;
    grp3.std = std(g3,'omitnan') ;
    
    % Put names in variable too
    grp1.names = gps.grp1_models ;
    grp2.names = gps.grp2_models ;
    grp3.names = gps.grp3_models ;

    % save
    sp = strcat(save_path,'AMOC Grouped Model Timeseries.mat') ;
    save(sp,'grp1','grp2','grp3','dt')

end % function