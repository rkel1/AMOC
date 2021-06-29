function [output,sWH,aWH] = function_apply_henderson(input,n_term) 

%%% Function to apply an n-term Henderson filter to an input timeseries.
%%% Based on code published by MATLAB and by Mark the Graph. Second and
%%% third output arguments give the symmetric and asymmetric weights
%%%
%%% Richard Kelson
%%% June 2020



L = length(input) ; % length of input timeseries 

if mod(n_term,1)
    n_term = ceil(n_term) ;
end
if ~mod(n_term,2) 
    n_term = n_term - 1 ;
end

% this function is coded below
[sWH,aWH] = function_henderson_weights(n_term) ;

% first and last terms
first = 1:n_term-1 ;
last = L-n_term+2:L ;

% apply filter
output = conv(input,sWH,'same') ;
output(L-floor(n_term/2)+1:end) = conv2(input(last),1,aWH,'valid');
output(1:floor(n_term/2)) = conv2(input(first),1,rot90(aWH,2),'valid');

end

function [sw,aw] = function_henderson_weights(n)

%%% Function to create the symmetrical [sw] and asymmetrical [aw]
%%% weightings for a Henderson filter of length n, where n is an odd
%%% number.
%%%
%%% Richard Kelson
%%% June 2020


%% Symmetrical
m  = (n-1)/2 ;
m1 = (m+1)^2 ;
m2 = (m+2)^2 ;
m3 = (m+3)^2 ;


% create the constant denominator
denominator = 8*sqrt(m2)*(m2-1)*(4*m2-1)*(4*m2-9)*(4*m2-25) ;

numerator = NaN(n,1) ;
% get the numerator
for i = 1:n
    i2 = (i-m-1)^2 ;
    numerator(i) = 315*(m1-i2)*(m2-i2)*(m3-i2)*(3*m2-11*i2-16) ;
end

sw = numerator./denominator ;

aw = zeros(n-1,m) ;
for i = 1:m
    asym_weight = henderson_asymmetric_weights(m+i,sw) ;
    aw(n-length(asym_weight):n-1,i) = asym_weight ;
end

aw = rot90(aw,2) ;

end

function aw = henderson_asymmetric_weights(m,w)

%%% Function to calculate asymmetric end weights [aw] for an input array
%%% of symmetrical weights [w] and the number of asymmetric weights [m].
%%%
%%% Richard Kelson
%%% June 2020

l = length(w) ;
assert(m < l,'[m] must be less than the length of [w]!')

% the beta squared / sigma squared 
ic = 1.0 ;
if l >= 13 && l < 15 
    ic = 3.5 ;
elseif l >= 15
 ic = 4.5 ;
end
b2s2 = (4/pi)/(ic^2) ;

% second term in formula
second = 0 ;
for j = m+1:l
    second = second + w(j)/m ;        
end

% third part (summation part)
third = 0 ;
for j = m+1:l
   third = third + (j - (m+1)/2)*w(j) ;
end

% denominator of third part
denom = 1 + m*(m-1)*(m+1)*b2s2 / 12 ;
    
aw = zeros(m,1) ;
for r = 1:m
    % numerator of third part
    numer = (r - (m+1)/2) * b2s2 ;

    % put all together
    aw(r) = w(r) + second + numer/denom * third ;
end

end
