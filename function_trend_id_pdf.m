function [segments,x,distrib] = ...
            function_trend_id_pdf(t,dt,timeseries,segment_in_months,x,...
                                  FILT_FLAG,n_term)

%%% Function to decompose the data timeseries into n segments of k months'
%%% length, and find the trend in those segments according to and OLS
%%% linear fit. Also fits a nonparametric, smoothed curve to the
%%% probability distribution function of trends, using a kernel density
%%% estimator such that the area under the curve is 1.
%%%
%%% The 5th input [x] gives the vector of x-values over which the fit
%%% should be calculated - if not given, this vector will default to
%%% equally-spaced points between the minimum and maximum of the trends
%%% found in the segments.
%%%
%%% The 6th input [FILT_FLAG] should only be changed from the default 0 if
%%% filtering should be performed AFTER segmentation has occurred. If a
%%% filter has a specified number of terms, this is picked in the 7th
%%% argument [n_term].
%%% Different types of filtering are performed based on the flag number:
%%%     1   -   n-term Henderson filtering
%%%     2   -   X-12 filtering
%%%
%%%
%%% Richard Kelson
%%% June 2020

% preallocate vector of trends
segments.trends = NaN(length(t)-segment_in_months,1) ;
segments.times = NaN(length(t)-segment_in_months,1) ;

if nargin < 6
    FILT_FLAG = 0 ;
    
    if nargin < 7
        n_term = 23 ;
    end
end

% get overlapping segments of k months
for i = 1:length(t)-segment_in_months 
    this_segment = timeseries(i:i+segment_in_months-1) ;
    this_dt      = dt(i:i+segment_in_months-1) ;   
    
    % get the time of the middle of the segment
    avg_time = mean(t(i:i+segment_in_months-1)) ;
        
    % filter this segment if requested
    switch FILT_FLAG
        case 0 % continue with unfiltered timeseries
            this_segment = this_segment ;
            
        case 1 % Henderson n-term filtering
            this_segment = function_apply_henderson(this_segment,n_term) ;
            
        case 2 % X-12 filtering
            this_segment = function_x12_filter(this_dt,this_segment) ;
    end
    
    % calculate the (linear) trend per year
    this_trend = fast_polyfit((1:segment_in_months)./12,this_segment,1) ;
    
    % put values into the struct
    segments.trends(i) = this_trend(1) ;
    segments.times(i)  = avg_time(1) ;
end

% detect if x values input, and if not, assigns default x
if nargin < 5
    x = linspace(min(segments.trends),max(segments.trends),20) ;
end

% fit a nonparametric, smoothed curve using a kernel density estimator with
% supersampling of x to get smoother fit
if nargout > 1
    distrib = pdf(fitdist(segments.trends,'kernel'),interp(x,5)) ;
    distrib = downsample(distrib,5) ;
end

% % return the supersampled x
% x = interp(x,5) ;

end