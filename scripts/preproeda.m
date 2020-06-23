%% PRE-PROCESS EDA
% Ivo Stuldreher
% modified on: 29/03/2019

% The function preproeda uses the Ledalab toolbox (http://www.ledalab.de/)
% to separate the phasic and tonic components in eda. Input argument (eda)
% is a struct that should at least contain fields 'data' and 'samplerate'.

function [phasic] = preproeda(conductance, samplerate)

%% Specify Ledalab format
% transfer eda structure in the required format for Ledalab processing
data = struct;
data.time = (1 : length(conductance))' / samplerate;
data.conductance = conductance;
data.timeoff = 0;

% create a struct event containing event information for Ledalab processing
time = [0]';
nid = [1]';
name = [];
userdata = [];
data.event = struct('time', time, 'nid', nid, 'name', name, 'userdata', userdata);

% add current path, save temporary matlab struct for Ledalab processing
scriptpath = pwd;
save(['Ledalab.mat'], 'data');

% process data using Ledalab, then return to current directory and load
% Ledalab output file. Seperate tonic and phasic components using
% continuous decomposition analysis (CDA)
date = datestr(now, 'yyyymmdd_hhMM');
Ledalab([pwd, '\Ledalab.mat'], ...
    'open', 'mat', ...
    'analyze', 'CDA', ...
    'overview', 0);
cd(scriptpath)
load(['Ledalab.mat'], 'analysis');
delete(['Ledalab.mat']);

% save phasic data in variable phasic
phasic = analysis.phasicData;


end