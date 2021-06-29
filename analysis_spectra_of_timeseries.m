% Script to create and display the spectra for CMIP6 models, as well as the
% RAPID array and the CMIP5 model mean.
%
% Richard Kelson
% June 2021

clear all

% file path parameters
home = pwd ;

cmip5_data_path = strcat(home,'/CMIP5/') ;
cmip6_data_path = strcat(home,'/AMOC 26 Data/') ;

save_data_path = sprintf('%s/Spectral Analysis/',home) ;
save_figure_path = sprintf('%s/Figures/Analysis/',home) ;

if ~exist(save_data_path,'dir')
    mkdir(save_data_path)
end

if ~exist(save_figure_path,'dir')
    mkdir(save_figure_path)
end



%% CMIP6 Analysis
% filtered
save_path = strcat(save_data_path,'CMIP6/Filtered/') ;
create_spectra_mat(cmip6_data_path,'%sAMOC_26_piControl_%s.mat',save_path,1)

% unfiltered
save_path = strcat(save_data_path,'CMIP6/Unfiltered/') ;
create_spectra_mat(cmip6_data_path,'%sAMOC_26_piControl_%s.mat',save_path,0)

%% CMIP5 Analysis    
% go through models, get spectra
save_path = strcat(save_data_path,'CMIP5/Filtered/') ;
create_spectra_mat(cmip5_data_path,'%s%s.mat',save_path,1)

% go through models, get spectra
save_path = strcat(save_data_path,'CMIP5/Unfiltered/') ;
create_spectra_mat(cmip5_data_path,'%s%s.mat',save_path,0)
%% Rapid Analysis
load('RAPID_trend.mat','times','timeseries','r_dt','var_inter') ;

% interpolate the raw timeseries to monthly data
r_amoc = interp1(times,timeseries,r_dt) ;

[rP_unf,rf_unf] = function_get_fft_spectra(r_amoc) ;
[rP_fil,rf_fil] = function_get_fft_spectra(var_inter) ;


%% Caesar analysis?


%% Plotting
fig = figure('Position',[100 471 1019 420]) ;
plt = tiledlayout(1,2) ;
title(plt,'Spectra of AMOC timeseries')

%%%% Unfiltered %%%%
nexttile(1) ; hold on
set(gca,'XScale','log','YScale','log','XDir','reverse')
grid on

ylim([1e-4 1e5])

xlabel('Period [years]')
ylabel('Power')
title('Unfiltered')

%%% plot CMIP6 models
c6_path = strcat(save_data_path,'CMIP6/Unfiltered/') ;
c6_files = dir(fullfile(c6_path,'*00*.mat')) ;
load(strcat(c6_path,c6_files.name))

% individual models
loglog(1./model_f(1,:),model_P,...
            'HandleVisibility','off',...
            'Color',[.6 .6 .6]) ;
        
% mean
loglog(1./model_f(1,:),mean(model_P,1,'omitnan'),...
            'DisplayName','CMIP6 mean',...
            'Color','b') ;


%%% plot CMIP5 model mean
c5_path = strcat(save_data_path,'CMIP5/Unfiltered/') ;
c5_files = dir(fullfile(c5_path,'*00*.mat')) ;
load(strcat(c5_path,c5_files.name))
        
% mean
loglog(1./model_f(1,:),mean(model_P,1,'omitnan'),...
            'DisplayName','CMIP5 mean',...
            'Color','r') ;

%%% plot RAPID
loglog(1./rf_unf,rP_unf,'DisplayName','RAPID','Color','k','LineWidth',1.25)



%%%% Filtered %%%%
nexttile(2) ; hold on
set(gca,'XScale','log','YScale','log','XDir','reverse')
grid on

ylim([1e-4 1e5])
xlim([1e0 2e1])

xlabel('Period [years]')
ylabel('Power')
title('Filtered, zoomed to periods of 1-20 years')
legend('location','eastoutside','box','off')

%%% plot CMIP6 models
c6_path = strcat(save_data_path,'CMIP6/Filtered/') ;
c6_files = dir(fullfile(c6_path,'*00*.mat')) ;
load(strcat(c6_path,c6_files.name))

% individual models
loglog(1./model_f(1,:),model_P,...
            'HandleVisibility','off',...
            'Color',[.6 .6 .6]) ;
        
% mean
loglog(1./model_f(1,:),mean(model_P,1,'omitnan'),...
            'DisplayName','CMIP6 mean',...
            'Color','b') ;


%%% plot CMIP5 model mean
c5_path = strcat(save_data_path,'CMIP5/Filtered/') ;
c5_files = dir(fullfile(c5_path,'*00*.mat')) ;
load(strcat(c5_path,c5_files.name))
        
% mean
loglog(1./model_f(1,:),mean(model_P,1,'omitnan'),...
            'DisplayName','CMIP5 mean',...
            'Color','r') ;

%%% plot RAPID
loglog(1./rf_fil,rP_fil,'DisplayName','RAPID','Color','k','LineWidth',1.25)

%% Save
sav_str = strcat(save_figure_path,'spectral_analysis.png') ;
exportgraphics(fig,sav_str) ;

%% Functions used in script
function create_spectra_mat(data_folder,search_string,save_path,FLAG)
    
    % flag is for deciding whether to filter
    
    if ~exist(save_path,'dir')
        mkdir(save_path)
    end

    % independent variable
    T = 12*498 ; % limit to first 498 years for comparison
    
    % see if we need to spceify the scenario
    if contains(search_string,'piControl') 
        models = get_model_names(data_folder,'piControl') ;
    else
        models = get_model_names(data_folder) ;
    end
    
    M = length(models) ;

    % preallocate for all_model variable
    model_f = NaN(M,T/2+1) ;
    model_P = NaN(M,T/2+1) ;

    for m = 1:M
        name = models{m} ;

        filename = sprintf(search_string,data_folder,name) ;

        load(filename,'dt','amoc') ;
        
        mT = length(dt) ; % length of loaded timeseries
        
        % restrict to maximum of 498 years
        if T < mT
            dt = dt(1:T) ;
            amoc = amoc(1:T) ;
        end

        % filter
        if FLAG
            amoc = function_x12_filter(dt,amoc) ;
        end

        % get spectrum
        [P,f] = function_get_fft_spectra(amoc) ;

        model_f(m,1:min(mT,T)/2+1) = f ;
        model_P(m,1:min(mT,T)/2+1) = P ;

        % save
        sav_str = strcat(save_path,name,'.mat') ;
        save(sav_str,'P','f','name') ;

    end

    % save all_model array
    sav_str = strcat(save_path,'00_ALL_MODEL.mat') ;
    save(sav_str,'model_P','model_f') ;

end
