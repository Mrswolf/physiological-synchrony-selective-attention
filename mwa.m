function [rho, rho_tt] = mwa(data1, data2, sr, varargin)
% MWA quantifies physiological synchrony (PS) using a moving window approach
%   [RHO, RHO_TT] = MWA(DATA1, DATA2) quantifies the similarity between
%   DATA1 and DATA2 using a moving window in which the Pearson correlation
%   is assessed. The moment-to-moment similarity is saved to RHO_TT, with
%   with an overall session index saved to RHO.
%
%   [...] = MWA(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies 
%   additional parameters and their values.  Valid parameters are the 
%   following:
%
%        Parameter         Value
%         'BartlettWindow'  1.5 (default) Bartlett running window size [s]
%         'CorWindow'   	15 (default) Pearson correlation running window
%                           size [s].
%         'CorStep'         1 (default) Pearson correlation running window
%                           step increment [s].
%         'SlopeWindow'     5 (default) moving slope running window size
%                           [s].
%         'SlopeStep'       1 (default) moving slope running window step
%                           increment [s].
%         'filtFlag'        false (default) / true filter data using 
%                           bartlett window
%         'slopeFlag'       false (default) / true differentiate data 
%                           before correlation calculation
%
%   See also: MOVSLOPE, MOVCORR
%
%   Copyright: Carl D. Marci, adopted by Ivo V. Stuldreher
%
%   Last modified on 15-04-2019
%       03-07-2019: Add flags for filter and slope, set standard to not
%       filter and differentiate the data.


%   This approach is obtained from [1]

%   REFERENCES
%   [1] Marci, C. D., Ham, J., Moran, E., & Orr, S. P. (2007). Physiologic 
%   correlates of perceived therapist empathy and social-emotional process 
%   during psychotherapy. The Journal of nervous and mental disease, 
%   195(2), 103-111.

% specify name value input arguments
p = inputParser;

defaultValue = 1.5;
errorMsg = 'Value must be positive, scalar and numeric.'; 
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), errorMsg);
addParameter(p, 'BartlettWindow', defaultValue, validationFcn);

defaultValue = 15;
errorMsg = 'Value must be positive, scalar and numeric.'; 
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), errorMsg);
addParameter(p, 'CorWindow', defaultValue, validationFcn);

defaultValue = 1;
errorMsg = 'Value must be positive, scalar and numeric.'; 
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), errorMsg);
addParameter(p, 'CorStep', defaultValue, validationFcn);

defaultValue = 5;
errorMsg = 'Value must be positive, scalar and numeric.'; 
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), errorMsg);
addParameter(p, 'SlopeWindow', defaultValue, validationFcn);

defaultValue = 1;
errorMsg = 'Value must be positive, scalar and numeric.'; 
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), errorMsg);
addParameter(p, 'SlopeStep', defaultValue, validationFcn);

defaultValue = false;
errorMsg = 'Value must be logical.'; 
validationFcn = @(x) assert(islogical(x), errorMsg);
addParameter(p, 'filtFlag', defaultValue, validationFcn);

defaultValue = false;
errorMsg = 'Value must be logical.'; 
validationFcn = @(x) assert(islogical(x), errorMsg);
addParameter(p, 'slopeFlag', defaultValue, validationFcn);

parse(p, varargin{:});
input.samplerate = sr;
input.data = [data1, data2];
input.time = (1 : length(data1))/sr;
    
% smoothen input data using a Bartlett window filter of length wt with step
% size st
if p.Results.filtFlag
    wl = p.Results.BartlettWindow * input.samplerate; % window length [num. samples]
    window = bartlett(wl);
    input.data = filter(window, 1, input.data);
end

% moving window slope calculation
if p.Results.slopeFlag
    input = movingslope(input, p.Results.SlopeWindow, p.Results.SlopeStep);
end

% moving window correlation calculation
rho = movingcorrelation(input, p.Results.CorWindow, p.Results.CorStep);
rho_tt = rho.data;
rho = log(sum(rho.data(rho.data > 0)) / sum(abs(rho.data(rho.data < 0))));

end

% ---------------------------------------------------------------------- %

function [slope] = movingslope(input, wt, st)

ws = wt * input.samplerate;
ss = st * input.samplerate;

% build filter coefficient to estimate the slope
p = mod(ws, 2);
s = (ws - p) / 2;
t = ((-s+1-p):s)';
slope.samplerate = st;

% calculate moving window slope with window size 'ws' and step increment of
% 'ss' by taking the mean of the gradient within 'ws'
slope.data = zeros(ceil(length(input.data)/ss), size(input.data,2));
slope.time = linspace(0, ceil(input.time(end)), length(slope.data))';
for i = 1 : size(input.data,2)
    for j = s+1 : ss : length(input.data) - s
        slope.data(ceil(j/ss), i) = mean(gradient(input.data(j-s:j+s, i)));
    end
    
end

end

function [rho] = movingcorrelation(input, wt, st)

ws = wt * input.samplerate;
ss = st * input.samplerate;

% build filter coefficient to estimate the slope
p = mod(ws, 2);
s = (ws - p) / 2;
t = ((-s+1-p):s)';
rho.samplerate = st;

% calculate moving window correlation with window size 'ws' and step increment of
% 'ss' by taking the mean of the gradient within 'ws'
rho.data = zeros(ceil(length(input.data)/ss), 1);
rho.time = linspace(0, ceil(input.time(end)), length(rho.data))';
for j = s+1 : ss : length(input.data) - s
    rho.data(ceil(j/ss)) = corr(input.data(j-s:j+s,1), input.data(j-s:j+s,2));
end

end
