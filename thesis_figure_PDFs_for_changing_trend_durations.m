% Script to plot 2 panels, showing how the PDF of trends changes when
% looking at trends of varying lengths. The first panel indicates how the
% observed trend at RAPID becomes more significant the longer it lasts; the
% second panel shows how the significance of the observed trend has changed
% in the last 5 years.
%
% Richard Kelson
% December 2020


clear all

% set flags
DRAFT_FLAG = 0 ; % add draft watermark to figures

% set file parameters
home = pwd ;
save_fig_path = sprintf('%s/Figures/Thesis Main/',home) ;
data_path = sprintf('%s/AMOC 26 Data/',home) ;

% set script parameters
scen = 'piControl' ;

% get model names from model_path directory
all_models = get_model_names(data_path,scen) ; 

%%% debug: don't use E3SM suite %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e3sm_idcs = contains(all_models,'E3SM') ;
models = all_models(~e3sm_idcs) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load RAPID data
load('RAPID_trend.mat','trend_var','tnd_yrs_var','yrs','duration') ;

T = 498*12 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

duration_effect_fig = figure('Position',[180 500 1250 460]) ;
plt = tiledlayout(1,2) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Panel 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load mean of model data
durs = round(12.*[10 duration 20 40]) ;  

edges = linspace(-0.3, 0.3, 300) ;

E = length(edges) ;
M = length(models) ;
D = length(durs) ;

% preallocate
freq       = NaN(E-1,M,D) ;
mu_array   = NaN(M,D) ;
std_array  = NaN(M,D) ;
dtb_array  = NaN(M,E,D) ;

for m_idx = 1:M
    
    model_name = models{m_idx} ;
    load_path = sprintf('%sAMOC_26_%s_%s.mat',data_path,scen,model_name) ;
    
    % check if mat file exists - if not, try to make it. If still bad, move
    % to next model
    if ~function_verify_mat_file(data_path,scen,model_name)
        continue
    end

    load(load_path,'amoc','dt','t') ;
    
    % take first 498 years only
    dt   = dt(1:T) ;
    t    = t(1:T) ;
    amoc = amoc(1:T) ;
    
    % get interannual component for varying annual cycle
    interannual = function_x12_filter(dt,amoc) ;
    
    % get PDF of given segment lengths
    for k_idx = 1:D 
        k = durs(k_idx) ;
        [pa,x,dtb] = function_trend_id_pdf(t,dt,interannual,k,edges) ;
        
        % find mean and standard deviation
        mu_array(m_idx,k_idx)  = mean(pa.trends) ;
        std_array(m_idx,k_idx) = std(pa.trends) ;
        
        % get frequencies
        n = histcounts(pa.trends,edges) ;
        
        % put in arrays
        freq(:,m_idx,k_idx) = n ;
        dtb_array(m_idx,:,k_idx) = dtb ;
    end
    
end
    
% get median/mean distribution for each segment length, and mean, std
median_dist = squeeze(median(freq,2,'omitnan')) ;
median_dtb  = squeeze(mean(dtb_array,1,'omitnan')) ;
median_mu   = median(mu_array,1,'omitnan') ;
median_std  = mean(std_array,1,'omitnan') ;

% get RAPID significance at each segment length
rapid_sig = abs(trend_var./median_std) ;

%% Write to .xls

%% Plot PDFs on one tile
% set colours
emph_c = [.7 .4 .5] ; 
demph_c = .6*[1 1 1] ;

cmap = NaN(D,3) ;
for d = 1:D
    if durs(d)==round(duration*12)
        cmap(d,:) = emph_c ;
    else
        cmap(d,:) = demph_c ;
    end
end

ax1 = nexttile(1) ; hold on


% plot de-emphasised durations
plot(x,median_dtb(:,durs~=round(duration*12)),'color',demph_c,'LineWidth',1) ; 

% plot RAPID duration 
plot(x,median_dtb(:,durs==round(duration*12)),'color',emph_c,'LineWidth',1) ; 

% plot the observed trend
plot(trend_var.*[1 1], [0 max(median_dtb(:))],'k-','LineWidth',1) ;
text(trend_var, .85*max(median_dtb(:)),{'Observed','trend after', '14.3 years'},...
        'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1])

% put on 2 std dev points for each 
pos = NaN(D,1) ;
for d = 1:D
    pos(d) = interp1(x,median_dtb(:,d),-2.*median_std(d)) ;
end
scatter(-2.*median_std,pos,25,cmap,'filled')

% add text labels to lines
str = num2str(round(durs'./12,1),'%.1f') ;
pos = interp1(x,median_dtb,0) ;

for d = 1:D
    t=text(0,pos(d),str(d,:),...
        'Color',cmap(d,:),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'BackgroundColor',[1 1 1],...
        'Margin',0.1) ;
end


% label axes
xlabel('AMOC Trend Value [Sv yr^{-1}]')
ylabel('Probability Density')
title('Significance if the observed trend persists')

xlim([-.3 .3])

% set panel lettering
xx = xlim() ;
yy = ylim() ;
text(xx(1) + 0.05.*abs(xx(1)),yy(2),'a)',...
        'FontSize',20,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top')
    
if DRAFT_FLAG % add watermark
    
    ax1.Units = 'pixels' ;
    ang = atand(ax1.Position(4)/ax1.Position(3)) ;
    
    t=text(0.5,0.65,'DRAFT',...
        'FontSize',60,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'Color',[.92 .92 .92],...
        'Rotation',ang,...
        'Units','Normalized') ;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Panel 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% get model mean

durs = round(12.*[10,11,12,13,duration]) ;  
edges = linspace(-0.45, 0.45, 300) ;

E = length(edges) ;
M = length(models) ; 
D = length(durs) ;

% preallocate
freq       = NaN(E-1,M,D) ;
mu_array   = NaN(M,D) ;
std_array  = NaN(M,D) ;
dtb_array  = NaN(M,5*E,D) ;

for m_idx = 1:M
    
    model_name = models{m_idx} ;
    load_path = sprintf('%sAMOC_26_%s_%s.mat',data_path,scen,model_name) ;
    load(load_path,'amoc','dt','t') ;
    
    % get interannual component for varying annual cycle
    interannual = function_x12_filter(dt,amoc) ;
    
    % get PDF of given segment lengths
    for k_idx = 1:D
        k = durs(k_idx) ;
        [pa,x,dtb] = function_trend_id_pdf(t,interannual,k,edges) ;
        
        % find mean and standard deviation
        mu_array(m_idx,k_idx)  = mean(pa.trends) ;
        std_array(m_idx,k_idx) = std(pa.trends) ;
        dtb_array(m_idx,:,k_idx) = dtb ;
        
        % get frequencies
        n = histcounts(pa.trends,edges) ;
        
        % put in frequency array
        freq(:,m_idx,k_idx) = n ;
    end
    
end
    
% get model-mean distribution for each segment length, and mean, std
median_dist = squeeze(median(freq,2,'omitnan')) ;
median_dtb  = squeeze(mean(dtb_array,1,'omitnan')) ;
median_mu   = median(mu_array,1,'omitnan') ;
median_std  = mean(std_array,1,'omitnan') ;

% get RAPID significance at each segment length
rapid_trends = interp1(round(yrs,3),tnd_yrs_var,round(durs/12,3)) ;
rapid_sig = abs(rapid_trends./median_std) ;

%% Plot on 1 tile
% set colours
emph_c = [.7 .4 .5] ; 
demph_c = .6*[1 1 1] ;

cmap = NaN(D,3) ;
for d = 1:D
    if durs(d)==round(duration*12)
        cmap(d,:) = emph_c ;
    else
        cmap(d,:) = demph_c ;
    end
end


nexttile(2) ; hold on

% plot PDFs
plot(x,median_dtb,'color',demph_c,'LineWidth',1) ; 
ax = gca ;

if DRAFT_FLAG % add watermark
    ax.Units = 'pixels' ;
    ang = atand(ax.Position(4)/ax.Position(3)) ;
    
    t=text(0.5,0.4,'DRAFT',...
        'FontSize',60,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'Color',[.92 .92 .92],...
        'Rotation',ang,...
        'Units','Normalized') ;
end

% find which is more negative: 2std dev or intercept? Plot lines between
% the two
for d = 1:D
    if -2*median_std(d) < rapid_trends(d) % not significant
        stt_val = -2.*median_std(d) ;
        end_val = rapid_trends(d) ;
        clr     = [0 0 .7] ;
    else % significant
        stt_val = rapid_trends(d) ;
        end_val = -2.*median_std(d) ;
        clr     = [.7 0 0] ;
    end
        
    
    x_seg = linspace(stt_val,end_val,50) ;
    segment = interp1(x,median_dtb(:,d),x_seg) ;
              
   plot(x_seg,segment,'Color',clr,'LineWidth',1)
end


% put on 2 std dev points for each 
s_pos = NaN(D,1) ;
for d = 1:D
    s_pos(d) = interp1(x,median_dtb(:,d),-2.*median_std(d)) ;
end
scatter(-2.*median_std,s_pos,20,demph_c,'filled')

% plot where the observations cut the PDF
o_pos = NaN(D,1) ;
for d = 1:D
    o_pos(d) = interp1(x,median_dtb(:,d),rapid_trends(d)) ;
end
scatter(rapid_trends,o_pos,20,'filled','k')


% add text labels to lines
str = num2str(round(durs'./12,1),'%.1f') ;
% str = {'2014' ; ...
%        '2015' ; 
%        '2016' ; ...
%        '2017' ; ...
%        '2018'} ;
   
pos = interp1(x,median_dtb,0) ;

for d = 1:D
    t = text(0,pos(d),str(d,:),...
        'Color',demph_c,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'BackgroundColor',[1 1 1],...
        'Margin',0.01,...
        'FontSize',10) ;
end

% set axes limits
xlim([min(x),max(x)])

% label axes
xlabel('AMOC Trend Value [Sv yr^{-1}]')
ylabel('Probability Density')
title('Significance if observations ended in 2014 to 2018')

% set panel lettering
xx = xlim() ;
yy = ylim() ;
text(xx(1) + 0.05.*abs(xx(1)),yy(2),'b)',...
        'FontSize',20,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top')
    

% create inset
ax2 = axes('Position',[0.7200 0.1500 0.1800 0.4500]) ;
box on 

idcs = find(x<-0.135 & x > -0.43 ) ; % expand this area of the graph

plot(x(idcs),median_dtb(idcs,:),'color',demph_c,'LineWidth',1) ; 
hold on


% put in text to show which years are being looked at
str1 = '2004 - 2014' ;
str2 = '2004 - 2018' ;

arw_x = [ax2.Position(1) + 0.5*ax2.Position(3), ...
         ax2.Position(1) + 0.8*ax2.Position(3)] ;
arw_y = [ax2.Position(3) + 0.35*ax2.Position(4), ...
         ax2.Position(3) + 0.1*ax2.Position(4)] ;

t1 = annotation('textarrow',arw_x, arw_y, ...
                    'String',str1,...
                    'units','normalized',...
                    'color',[.2 .2 .2]) ;
                
t2 = text(ax2,.99,.05,str2,'units','normalized',...
                           'HorizontalAlignment','right',...
                           'VerticalAlignment','bottom',...
                           'color',[.2 .2 .2]) ;

% put on 2 std dev points for each 
scatter(-2.*median_std,s_pos,20,demph_c,'filled')

% plot where the observations cut the PDF
scatter(rapid_trends,o_pos,20,'filled','k')

% plot lines between 2std dev point and obs intercept
for d = 1:D
    if -2*median_std(d) < rapid_trends(d) % not significant
        stt_val = -2.*median_std(d) ;
        end_val = rapid_trends(d) ;
        clr     = [0 0 .7] ;
    else % significant
        stt_val = rapid_trends(d) ;
        end_val = -2.*median_std(d) ;
        clr     = [.7 0 0] ;
    end
        
    
    x_seg = linspace(stt_val,end_val,50) ;
    segment = interp1(x,median_dtb(:,d),x_seg) ;
              
   plot(x_seg,segment,'Color',clr,'LineWidth',1)
end
  

                  
% set inset axis formatting
xlim([min(x(idcs)),max(x(idcs))]) 
ylim([.05 1.1*max([o_pos; s_pos])])
set(ax2,'xtick',[],'xticklabels',[])
set(ax2,'ytick',[],'yticklabels',[])

% get inset limits for indication on main plot
ax2_x = xlim ;
ax2_y = ylim ;
rect_pos = [ax2_x(1), ax2_y(1), diff(ax2_x), diff(ax2_y)] ;

% return to main axis, put a box around magnified section
rectangle(ax,'Position',rect_pos)
plot(ax, [rect_pos(1) +    rect_pos(3),   0], ...
         [rect_pos(2) + .5*rect_pos(4), 1.5], ...
         'k--')
     



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DRAFT_FLAG
    sp = sprintf('%sDRAFT-PDFs affected by Duration.png',save_fig_path) ;
else
    sp = sprintf('%sPDFs affected by Duration.png',save_fig_path) ;
end
exportgraphics(duration_effect_fig,sp) ;