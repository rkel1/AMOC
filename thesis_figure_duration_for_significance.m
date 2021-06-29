% Script to plot a figure showing how a given trend gets more significant
% (i.e. reaches further from a distribution's mean) if it lasts for a
% longer amount of time; also plots the standard deviation of each model's
% PDF at different times. Done for each model.
%
% Richard Kelson
% January 2020

clear all

fprintf(1,'Running %s:\n%s\n',mfilename,repmat('=',1,70)) ;

% file parameters
home = pwd ;
load_data_path = sprintf('%s/AMOC Sig vs Length/',home) ;
save_fig_path = sprintf('%s/Figures/Thesis Main/',home) ;

if ~exist(save_fig_path,'dir')
    mkdir(save_fig_path)
end

% load RAPID data
load('RAPID_trend.mat','duration','trend_var','tnd_yrs_var','yrs') ;

% load the CMIP6 distributions files
mat_file = strcat(load_data_path,'VAC_model_VAC_rapid_sig_vs_length.mat') ;
load(mat_file,'trend_lengths_indep',...
              'std_dev_of_PDFs',...
              'sigma_array',...
              'lengths',...
              'models') ;
vac = load(mat_file) ;
          
% load the CAC CMIP6 distributions file for comparison
mat_file = strcat(load_data_path,'CAC_model_CAC_rapid_sig_vs_length.mat') ;
cac = load(mat_file) ;

% label representative models
r_mds = {'ACCESS-CM2','CNRM-ESM2-1','NorESM2-MM'} ;
m_idcs = find(contains(models,r_mds)) ;

%% get 'emergence time' for each models (duration to hit 2 std dev)
M = length(models) ;
emergence_time = NaN(M,1) ;

for m = 1:M
    emergence_time(m) = interp1(sigma_array(m,:),trend_lengths_indep,2) ;
end

% normalise
pivot = duration/max(abs(emergence_time)) ; % for getting white at 14.3 years
n_emg_t = emergence_time./max(abs(emergence_time)) ;

%% Fit std dev lines to decaying exponential
[fit_coeffs, resid_vals] = fit_model(vac) ;

%% Get varying-constant emergence time difference for each model
emergence_time_cte = NaN(M,1) ;

for m = 1:M
    emergence_time_cte(m) = interp1(vac.sigma_array(m,:),...
                                    vac.trend_lengths_indep,2) ;
end

% calculate the difference
emg_diffs = emergence_time_cte - emergence_time ;

%% get colourmap for each model (based on emergence time)
% load colourmap
cmap = cmocean('balance','pivot',pivot) ; close all
Y = linspace(0,1,256) ;

colours = NaN(M,3) ;
for m = 1:M
    colours(m,:) = interp2(1:3, Y, cmap, 1:3, n_emg_t(m)) ;
    
    % set colour if no emergence time
    if all(isnan(colours(m,:)))
        colours(m,:) = [0.2 0 0] ;
    end
end

save('cmap_emergence_time.mat','colours','models') ;

%% Setup figure
duration_fig = figure('Position',[200 150  1109 413]) ;
plt = tiledlayout(1,2) ;
% Plot 2x Standard Deviation change per duration per model
ax = nexttile(1) ;
hold on
set(gca,'XScale','log','YScale','log')


xlabel('Segment length [years]')
ylabel('2\sigma from model PDF [Sv yr^{-1}]')
title('Emergence time across models')

ylim([0.01 1.1])

% per-model
for m = 1:M
    plot(trend_lengths_indep,2*std_dev_of_PDFs(m,:),'Color',colours(m,:))
end

% mean
plot(trend_lengths_indep,mean(2*std_dev_of_PDFs,1),...
    '-',...
    'Color',.2*[1 1 1],...
    'LineWidth',2)

% emergence line - changing RAPID trend value
plot(yrs,abs(tnd_yrs_var),'k--','LineWidth',1.25)

% % text label for observation line
% str = '0.15 Sv yr^{-1}' ;
% text(5.4,abs(trend_var)+.006,str,...
%      'BackgroundColor','w',...
%      'VerticalAlignment','middle',...
%      'FontSize',10)

% labelling
% set panel lettering
text(.97,.99,'a)',...
        'FontSize',20,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top',...
        'Units','Normalized',...
        'Color','k')
    
 

% inset
% create inset
ax2 = axes('Position',[0.125 0.15 0.1500 0.450]) ;
box on 
hold on
set(gca,'XScale','log','YScale','log')

idcs = find(trend_lengths_indep>=15 & trend_lengths_indep < 18) ; % expand this area of the graph

for m = 1:M
    loglog(trend_lengths_indep(idcs),2*std_dev_of_PDFs(m,idcs),'Color',colours(m,:))
end

ylim([.075 .278])

% text labels (matching panel d)
xvals = linspace(15.25,16.75,5) ;
x_count = 1 ;

% sort by value at 10.2
[y_vals,s_idcs] = sort(interp2(1:M,...
                             trend_lengths_indep(idcs),...
                             2*std_dev_of_PDFs(:,idcs)',...
                             1:M,...
                             xvals(x_count)),...
                             'descend') ;

for m = 1:M
    % go through models in decreasing value of sigma
    s = s_idcs(m) ;
    s_m = models{s} ;
    
    yval = interp1(trend_lengths_indep(idcs),...
                   2*std_dev_of_PDFs(s,idcs),...
                   xvals(x_count)) ;
    
    % skip if outside limits
    if yval > max(ylim(ax2))
        continue
    end
    
    % emphasise if one of repr models
    if any(contains(r_mds,s_m))
        txt = strcat('\color[rgb]{0.6914,0.2617,0.2617}\bf',num2str(s)) ;
    else
        txt = num2str(s) ;
    end
    
    t=text(xvals(x_count),yval,txt,...
         'VerticalAlignment','middle',...
         'HorizontalAlignment','center',...
         'BackgroundColor',[1 1 1],...
         'FontSize',7,...
         'Margin',0.001) ;
               
    % increase x_count          
    x_count = x_count + 1 ;
    if x_count > length(xvals)
        x_count = 1 ;
    end
    
end

% get rid of X,Y ticks
set(ax2,'XTick',[],'YTick',[])

% get x and y limits of box
xlims = xlim(gca) ;
ylims = ylim(gca) ;

rect_pos = [xlims(1), ylims(1), diff(xlims), diff(ylims)] ;

% return to main axes, put a box around magnified section
set(duration_fig,'CurrentAxes',ax)
rectangle(ax,'Position',rect_pos)
plot(ax, [rect_pos(1) + 0.0*rect_pos(3),   11], ...
         [rect_pos(2) + 0.4*rect_pos(4), .035], ...
         'k-')
     
% Plot outlier on main axis
s = s_idcs(1) ;
s_m = models{s} ;
xval=22 ;
yval = interp1(trend_lengths_indep,...
                   2*std_dev_of_PDFs(s,:),...
                   xval) ;



t=text(ax,xval,yval,num2str(s),...
         'VerticalAlignment','middle',...
         'HorizontalAlignment','center',...
         'BackgroundColor',[1 1 1],...
         'FontSize',9,...
         'Margin',0.001) ;

%% Plot the fit coefficients and emergence time
nexttile(2) ;
hold on

% axis equal
ylim([0 .3])
xlim([0 .14])

slope = 1/(range(ylim)/range(xlim)) ;

% C
x = fit_coeffs(:,3).*ones(M,2) ;


alpha = abs(trend_var) ; % significance threshold

% A*tau^(-b)
tau = 14.3 ; % corresponds to 15 years
y1 = fit_coeffs(:,1).*tau.^(-fit_coeffs(:,2)) ;

tau = 14.3+10 ; % corresponds to 15 years
y2 = fit_coeffs(:,1).*tau.^(-fit_coeffs(:,2)) ;

y = horzcat(y1, y2) ;

% plot trend lines
rapid_trends = [0.3224, 0.2668, 0.2101, 0.1482] ;
R = length(rapid_trends) ;


for r = 1:R 
    clr = [.6 .6 .6] ;
    width = 0.5 ;
    
    if r == R
        clr = [0 0 0] ;
        width = 2 ;
    end
    
    plot([0 rapid_trends(r)], [rapid_trends(r) 0], ...
         'LineWidth',width,...
         'Color',clr)    
end



% plot a 45-degree line
endpoint = .19 ;
plot([0 endpoint], [0 endpoint],'k--')

% add text to line
pos = 0.117 ;
t = text(pos,pos+0.003,'At^{-b} = C',...
        'Color','k',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle' ,...
        'BackgroundColor','w',...
        'Margin', 2,...
        'FontSize',9) ;
set(t,'Rotation',45*slope) ;


% plot models
scatter(x(:,1),y1,25,colours,'filled') 
scatter(x(:,1),y2,5,colours,'filled') 
for m = 1:M
    plot(x(m,:)',y(m,:)','Color',colours(m,:),'LineWidth',1)
end

% add text to trend lines
xvals = [0.07, 0.06, 0.05, 0.04] ; % x values for text
for r = 1:R 
    clr = [.6 .6 .6] ;
    if r == R
        clr = [0 0 0] ;
    end
    
    str = sprintf('%.2f Sv yr^{-1}',rapid_trends(r)) ;
    
    xt = linspace(0, 0.3, 100) ;
    yt = linspace(0.3, 0, 100) ;
    yval = interp1(xt,yt,xvals(r)) ;
    
    t = text(xvals(r),rapid_trends(r)-xvals(r)+0.0025,str,...
        'Color',clr,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle' ,...
        'Margin', .6,...
        'BackgroundColor','w') ;
        set(t,'Rotation',-40*slope) ;
end

% labelling
xlabel('C')
ylabel('At^{-b}')
title(strcat('Short- vs. long-term variability: t = [14.3 24.3]'))

text(x(m_idcs(1),1)+0.001,y1(m_idcs(1))+0.004,'\leftarrow ACCESS-CM2',...
     'HorizontalAlignment','left',...
     'FontSize',8,...
     'BackgroundColor','w',...
     'Margin',0.1)

text(x(m_idcs(2),1)+0.001,y1(m_idcs(2)),'\leftarrow CNRM-ESM2-1',...
     'HorizontalAlignment','left',...
     'FontSize',8,...
     'BackgroundColor','w',...
     'Margin',0.1)

text(x(m_idcs(3),1)+0.002,y1(m_idcs(3))+0.004,'\leftarrow NorESM2-MM',...
     'HorizontalAlignment','left',...
     'FontSize',8,...
     'BackgroundColor','w',...
     'Margin',0.01)
 
% label outlier with asterisk
o_idx = find(contains(models,'FGOALS-g3')) ;
text(x(o_idx)+0.001,y2(o_idx)-0.002,'*',...
    'HorizontalAlignment','left',...
    'FontSize',14)

% set panel lettering
text(.97,.99,'b)',...
        'FontSize',20,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top',...
        'Units','Normalized',...
        'Color','k')

%% Plot the fit coefficients and emergence time
% nexttile(3) ;
% hold on
% 
% % load variables
% sig_pi_path = sprintf('%s/AMOC Sig vs Length/',home) ;
% sig_pd_path = sprintf('%s/AMOC Sig vs Length (Hist)/',home) ;
% met_path = sprintf('%s/AMOC Timeseries Variation/',home) ;
% 
% 
% sig_pi_struct = load(sprintf('%sVAC_model_VAC_rapid_sig_vs_length.mat',sig_pi_path)) ;
% sig_pd_struct = load(sprintf('%sVAC_model_VAC_rapid_sig_vs_length.mat',sig_pd_path)) ;
% load('RAPID_trend.mat','duration','trend_cte','var_inter') ; 
% 
% ma_pi = load('Metric Array PI.mat') ;
% ma_pd = load('Metric Array PD.mat') ;
% 
% % Get standard deviation of each model's PD PDF
% s14 = function_get_std_dev(models,sig_pd_struct,duration) ;
% 
% % Get model mean std dev of PI PDFs and monthly mean
% s14_pi_mean = mean(function_get_std_dev(models,sig_pd_struct,duration)) ;
% s14_pi_mm   = mean(ma_pi.metric_array(1:M,9)) ;
% 
% % Get correlation across models for each metric in PD
% [r14,p14] = function_get_corr(ma_pd.metric_array,s14) ;
% 
% % Plot individual models
% scatter(ma_pd.metric_array(1:M,9),s14,25,colours,'filled','HandleVisibility','off') ;
% 
% % Plot model means
% scatter(mean(ma_pd.metric_array(1:M,9)),mean(s14),25,'k','filled','DisplayName','CMIP6 PD Mean')
% scatter(s14_pi_mm,s14_pi_mean,25,'d','k','DisplayName','CMIP6 PI Mean')
% 
% % Plot Caesar value
% load('Caesar_AMOC.mat')
% std_caes = std(amoc) ;
% trends14 = function_trend_id_pdf(t,dt,amoc,round(duration*12)) ;
% std_caes_tnds = std(trends14.trends) ;
% scatter(std_caes,std_caes_tnds,25,'k*',...
%         'DisplayName','Caesar et al.')
%     
% % Plot CMIP5 value
% cmip5 = load(strcat(home,'/CMIP5 Analysed Data/00_ALL_MODEL_analysis.mat')) ;
% cmip5_14_year_point = mean(interp2(cmip5.trend_lengths_indep,1:17,...
%                                     cmip5.model_std_of_PDFs,...
%                                     round(duration*12)/12,1:17)) ;
% cmip5_pt = [mean(cmip5.model_std_amoc) cmip5_14_year_point] ;
% scatter(cmip5_pt(1),cmip5_pt(2),25,'k','DisplayName','CMIP5 PI mean') ;
% 
% 
% % Plot RAPID x-value (std dev of timeseries)
% plot(std(var_inter).*[1 1], [0 .3],...
%      'k--','LineWidth',2,'DisplayName','RAPID observations') ;
% 
% % label representative models
% text(ma_pd.metric_array(m_idcs(1),9)-0.02,s14(m_idcs(1))+.0015,'ACCESS-CM2 \rightarrow',...
%      'HorizontalAlignment','right',...
%      'VerticalAlignment','middle',...
%      'FontSize',8,...
%      'BackgroundColor','None',...
%      'Margin',0.1)
%  
% text(ma_pd.metric_array(m_idcs(2),9)+0.03,s14(m_idcs(2))+.007,'CNRM-ESM2-1',...
%      'HorizontalAlignment','left',...
%      'VerticalAlignment','middle',...
%      'FontSize',8,...
%      'BackgroundColor','None',...
%      'Margin',0.1)
%  
% text(ma_pd.metric_array(m_idcs(3),9)-0.03,s14(m_idcs(3))+.0015,'NorESM2-MM \rightarrow',...
%      'HorizontalAlignment','right',...
%      'VerticalAlignment','middle',...
%      'FontSize',8,...
%      'BackgroundColor','None',...
%      'Margin',0.1)
%  
% 
% % labelling
% ylabel(sprintf('\\sigma(%.1f-year AMOC trend) [Sv yr^{-1}]',duration))
% xlabel('\sigma(AMOC) [Sv]')
% title('Monthly vs. 14.3-year trend variability')
% 
% % formatting
% ylim([0 .3])
% xlim([.5 2.6])
% 
% % legend
% lg = legend('location','northwest','box','off') ;
% 
% % set panel lettering
% text(.01,.99,'c)',...
%         'FontSize',20,...
%         'HorizontalAlignment','left',...
%         'VerticalAlignment','top',...
%         'Units','Normalized',...
%         'Color','k')
    
%% Chart the effect of using the X-12 algorithm
% nexttile(4)
% hold on
% 
% b = barh(emg_diffs,'FaceColor','flat') ;
% b.CData = colours ;
% b.FaceAlpha = 0.8 ;
% b.EdgeAlpha = 0 ;
% 
% 
% 
% % label bars with models
% m_labels = cell(M,1);
% for m = 1:M
%     m_labels{m} = strcat('\color{black} ',sprintf('%.0f. %s',m,models{m})) ; 
%     if any(contains(r_mds,models{m}))
%         m_labels{m} = strcat('\color[rgb]{0.6914,0.2617,0.2617}\bf',sprintf('%.0f. %s',m,models{m})) ;
%     end 
% end
% 
% yticks(1:M) 
% yticklabels(m_labels) 
% ax = gca ;
% ax.FontSize = 8 ;
% set(gca,'YDir','reverse')
% 
% % remove unecessary lines
% set(gca,'XTick',[])
% b.Parent.XAxis.Color = [1 1 1] ;
% 
% % flip plot horizontally
% b.Parent.YAxisLocation = 'right' ;
% set(gca,'XDir','reverse')
% 
% % direct label bars
% for m = 1:M
%     clr = [0 0 0] ;
%     if norm(b.CData(m,:)) < .75
%         clr = [1 1 1] ;
%     end
%     
%     text(b.YData(m)-0.025, m, sprintf('%.1f',b.YData(m)),...
%         'FontSize',7,...
%         'Color',clr,...
%         'HorizontalAlignment','left')
% end
% 
% 
% % labelling
% title('Decrease in emergence time due to X-12 (years)','FontSize',11.5)
% 
% % set panel lettering
% t = text(.02,.99,'d)',...
%         'FontSize',20,...
%         'HorizontalAlignment','left',...
%         'VerticalAlignment','top',...
%         'Units','Normalized',...
%         'Color','k') ;

    
%% Save colourmap information
cmap = cell(M,2) ;
cmap(:,1) = models ;

for m = 1:M
    cmap(m,2) = {colours(m,:)} ;
end
save('model_cmap.mat','cmap')
    
%% Save figure
sp = strcat(save_fig_path,'Duration for Significance.png') ;
exportgraphics(duration_fig,sp) ;


%% Extra (optional) plot
% collapse_fig = figure() ;
% 
% ax = nexttile(1) ;
% hold on
% set(gca,'XScale','log','YScale','log')
% 
% % per-model
% for m = 1:M
%     plot(trend_lengths_indep,2*std_dev_of_PDFs(m,:) - fit_coeffs(m,3),'Color',colours(m,:))
% end
% 
% 
% xlabel('Segment length [years]')
% ylabel('2\sigma - C [Sv yr^{-1}]')
% title('2\sigma less C')
% 
% 
% ylim([0.01 1.1])
% 
% sp = strcat(home,'/Figures/Trend Analysis/','variation less C.png') ;
% exportgraphics(collapse_fig,sp) ;


%% Functions 
function std_devs = function_get_std_dev(models,sig_struct,duration)

    M = length(models) ;

    std_devs = NaN(M,1) ;

    for m = 1:M
        name = models{m} ;

        % get index in sig trend relating to this model
        s_idx = find(strcmp(sig_struct.models,name)) ;

        % if empty, skip
        if isempty(s_idx)
            continue
        end

        % get value of 1 std dev at 14.4 years, in Sv/yr
        sig = interp1(sig_struct.trend_lengths_indep,sig_struct.std_dev_of_PDFs(s_idx,:),duration) ;
        std_devs(m) = sig ;
    end

end

function [r,p] = function_get_corr(metric_array,std_devs,bad_idcs)
    
    T = size(metric_array,2) ;
    M = length(std_devs) ;
    
    if ~exist('bad_idcs','var')
        bad_idcs = zeros(M,1) ;
    end
    
    idcs = 1:M ;
    idcs = idcs(~bad_idcs) ;

    % Get correlation across models for each metric
    r = NaN(T,1) ;
    p = NaN(T,1) ;

    for ii = 1:T
        [rii,pii] = corrcoef(metric_array(idcs,ii),std_devs(idcs),'rows','complete') ;
        r(ii) = rii(2) ;
        p(ii) = pii(2) ;
    end
end
