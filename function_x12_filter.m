function [interannual,seasonal,subannual,dt,t] = ...
    function_x12_filter(times, timeseries, n_term)

%%% Script to apply the X-12 timeseries decomposition methodology detailed
%%% in MATLAB's Econometrics Toolbox Documentation (section detailing
%%% Seasonal Adjustment). Of use when the annual cycle of a timeseries
%%% varies. Input times must be a datetime object
%%%
%%% The input data timeseries is decomposed into interannual, seasonal and
%%% subannual components that have a monthly timestep. These monthly
%%% timesteps are also output as a datetime vector (dt) and a double vector
%%% (t) - the latter values are the year of the timestep.
%%%
%%% Optional arguments:
%%%     n_term - to hard code the number of terms to use in the Henderson
%%%                 filter
%%%
%%% Richard Kelson
%%% June 2020

if nargin < 3
    n_term = 23 ; % number of terms in henderson filter - must be odd
end
    
%% check inputs in a consistent length format
assert(length(times) == length(timeseries),...
    'Input times must be the same length as the input timeseries')

assert(isa(times,'datetime'),...
    'Input times must be given as a datetime object')

%% check whether input times are monthly
[monthly_series,dt,t] = make_monthly(timeseries,times) ;


%% Start algorithm
L = length(monthly_series) ; % length of timeseries

% remove trend using a moving average specified by cutoff
s = 12 ;
sW_mva = [1/(2*s); repmat(1/s,(s-1),1);1/(2*s)] ;
data_trend = conv(monthly_series,sW_mva,'same') ;

% to avoid problems at beginning and end of data with moving average, set
% these values to the closest value that is fully calculated
data_trend(1:s/2) = data_trend(s/2+1) ;
data_trend(L-s/2+1:L) = data_trend(L-s/2) ;

data_detrend_1 = monthly_series - data_trend ;

% Apply S(3,3) filter
data_seasonal_1 = function_apply_seasonal_filter(data_detrend_1,3,3) ;

% iterate: deseason series
data_deseason_1 = monthly_series - data_seasonal_1 ;

% apply h_term henderson filter
data_hend_trend_1 = function_apply_henderson(data_deseason_1,n_term) ;

% detrend timeseries
data_detrend_2 = monthly_series - data_hend_trend_1 ;

% Apply S(3,5) filter
data_seasonal_2 = function_apply_seasonal_filter(data_detrend_2,3,5) ;


% save to output variables
interannual = data_hend_trend_1 ;
seasonal = data_seasonal_2 ;
subannual = monthly_series - interannual - seasonal ;

end