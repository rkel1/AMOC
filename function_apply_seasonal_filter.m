function output = function_apply_seasonal_filter(input,n,m)

% Function to apply a S(n,m) seasonal filter to an input timeseries.
% Assumes calculating for an annual cycle where data has monthly timesteps.

L = length(input) ;

% store indices for each month (i.e. all the Jans, then Febs etc.)
s = 12 ;
s_idx = cell(s,1) ;
for i = 1:s
    s_idx{i,1} = i:s:L ;
end

% weights for 2x12 moving average
sW13 = [1/24; repmat(1/12,11,1);1/24] ;

% get symmetric weights
matA = repmat(1/n,n,1) ;
matB = repmat(1/m,m,1) ;
sw = conv(matA,matB) ;

% get asymmetric weights
switch m
    case 3
        aw = [.259 .407; 
              .37  .407;
              .259 .185;
              .111    0];
    case 5
        aw =  [.150 .250 .293;
               .217 .250 .283;
               .217 .250 .283;
               .217 .183 .150;
               .133 .067    0;
               .067   0     0];
    otherwise
        error('A S(%g,%g) filter is not supported!',n,m)
end

% apply filter
seas = NaN*input ;

for i = 1:s
    ns = length(s_idx{i}) ;
    first = 1:m+1 ;
    last = ns-m:ns ;
    dat = input(s_idx{i}) ;
    
    sd = conv(dat,sw,'same') ;
    sd(1:size(aw,2))  = conv2(dat(first),1,rot90(aw,2),'valid') ;
    sd(ns-(size(aw,2)-1):ns) = conv2(dat(last),1,aw,'valid') ;
    seas(s_idx{i}) = sd ;
end

% centre seasonal component around zero by finding trend
sb = conv(seas,sW13,'same') ;
sb(1:6) = sb(7) ;
sb(L-5:L) = sb(L-6) ;

output = seas - sb ;

end