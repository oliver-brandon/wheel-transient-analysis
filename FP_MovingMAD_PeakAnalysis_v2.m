clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = "2";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_MovingMAD_PeakAnalysis.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% input paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 20; % time threshold below which we will discard
N = 10; %Downsample N times
SAMPLE_RATE = 1017/N;
SESSION_DURATION = 3600;
WINDOW_SIZE_SECONDS = 15;
MAD_MULTIPLIER = 2;
MIN_PK_WIDTH = 0.2;
%Snippet Args%
SNIPPET_DURATION = 10; % Duration of the snippet in seconds
SNIPPET_START_TIME = 1000; % Start time of the snippet in seconds
%Stream Stores%
DLS_ISOS = 'x405A'; % name of the 405A store
DLS_DA = 'x465A'; % name of the 465A store
NAC_ISOS = 'x405C'; % name of the 405C store
NAC_DA = 'x465C'; % name of the 465C store

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Version: %s\n",VERSION)
myDir = uigetdir('H:\My Drive\wheel-peak-test\tanks');
BLOCKPATH = myDir;
data = TDTbin2mat(BLOCKPATH, 'T2', SESSION_DURATION, 'TYPE', {'streams'});

%time array used for all streams%
time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%

ind = find(time>T,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
data.streams.(NAC_DA).data = data.streams.(NAC_DA).data(ind:end);
data.streams.(NAC_ISOS).data = data.streams.(NAC_ISOS).data(ind:end);

%downsample streams and time array by N times%
data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
data.streams.(NAC_ISOS).data = downsample(data.streams.(NAC_ISOS).data, N);
data.streams.(NAC_DA).data = downsample(data.streams.(NAC_DA).data, N);
time = downsample(time, N);

%detrend & dFF%
%465A%
bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
detrend_465A = detrend(dFF);
detrend_465A = zscore(detrend_465A);
%465C%
bls2 = polyfit(data.streams.(NAC_ISOS).data,data.streams.(NAC_DA).data,1);
Y_fit_all2 = bls2(1) .* data.streams.(NAC_ISOS).data + bls2(2);
Y_dF_all2 = data.streams.(NAC_DA).data - Y_fit_all2; %dF (units mV) is not dFF
dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
std_dFF2 = std(double(dFF2));
detrend_465C = detrend(dFF2);
detrend_465C = zscore(detrend_465C);

%set up moving MAD window%

% determines if window_size evenly divides into the streams and pads with 
% NaN if not
window_size = ceil(WINDOW_SIZE_SECONDS * SAMPLE_RATE);
remainderA = ceil((window_size - mod(length(detrend_465A),window_size)));
remainderC = ceil((window_size - mod(length(detrend_465C),window_size)));
if remainderA == window_size
    signalA = detrend_465A;
    signalA_ind = 1:length(signalA);
else
    padding = NaN(1,remainderA);
    signalA = [detrend_465A,padding];
    signalA_ind = 1:length(signalA);
end

if remainderC == window_size
    signalC = detrend_465C;
    signalC_ind = 1:length(signalC);
else
    padding = NaN(1,remainderC);
    signalC = [detrend_465C,padding];
    signalC_ind = 1:length(signalC);
end
time_pad = (1:length(signalA))/SAMPLE_RATE;

signalA_pks = [];
for i = 1:window_size:length(signalA)
    if (i + window_size - 1) > length(signalA)
        window_end = length(signalA);
    else
        window_end = i + window_size - 1;
    end
    
    sigA_med = median(signalA(i:window_end));
    sigA_mad = mad(signalA(i:window_end), 1);
    sigA_thr = (sigA_med + (MAD_MULTIPLIER * sigA_mad));

    [pksA, locsA] = findpeaks(signalA(i:window_end), ...
        time_pad(i:window_end), 'MinPeakHeight', sigA_thr, ...
        'MinPeakDistance', MIN_PK_WIDTH);

    signalA_pks = [signalA_pks, locsA];
end

signalC_pks = [];
for i = 1:window_size:length(signalC)
    if (i + window_size - 1) > length(signalC)
        window_end = length(signalC);
    else
        window_end = i + window_size - 1;
    end
    
    sigC_med = median(signalC(i:window_end));
    sigC_mad = mad(signalC(i:window_end), 1);
    sigC_thr = (sigC_med + (MAD_MULTIPLIER * sigC_mad));
    % sigC_thr = madMultiplier * sigC_mad;

    [pksC, locsC] = findpeaks(signalC(i:window_end), ...
        time_pad(i:window_end), 'MinPeakHeight', sigC_thr, ...
        'MinPeakDistance', MIN_PK_WIDTH);

    signalC_pks = [signalC_pks, locsC];
end

x1 = ceil(time(1,1));
x2 = ceil(time(1,end));
% Plotting
f1 = figure;
plot(time_pad, signalA)
xlim([x1 x2]);
hold on
peak_indices = ismember(time_pad, signalA_pks);
plot(time_pad(peak_indices), signalA(peak_indices), 'ro')
hold off

f2 = figure;
plot(time_pad, signalC)
xlim([x1 x2]);
hold on
peak_indices = ismember(time_pad, signalC_pks);
plot(time_pad(peak_indices), signalC(peak_indices), 'ro')
hold off


%% Snippet Plotting %%
% Find the indices corresponding to the start and end times of the snippet
[~, snippet_start_index] = min(abs(time_pad - SNIPPET_START_TIME));
snippet_end_index = snippet_start_index + round(SNIPPET_DURATION / (time_pad(2) - time_pad(1))) - 1;

% Ensure the end index does not exceed the length of the signal
snippet_end_index = min(snippet_end_index, length(signalA));

% SIGNAL A %
% Plot the snippet of the signal
f3 = figure;
plot(time_pad(snippet_start_index:snippet_end_index), signalA(snippet_start_index:snippet_end_index))
hold on

% Plot the peaks within the snippet
peak_indices_within_snippet = peak_indices & (time_pad >= SNIPPET_START_TIME) & (time_pad <= (SNIPPET_START_TIME + SNIPPET_DURATION));
plot(time_pad(peak_indices_within_snippet), signalA(peak_indices_within_snippet), 'ro')
hold off

% SIGNAL C %
% Plot the snippet of the signal
f4 = figure;
plot(time_pad(snippet_start_index:snippet_end_index), signalC(snippet_start_index:snippet_end_index))
hold on

% Plot the peaks within the snippet
peak_indices_within_snippet = peak_indices & (time_pad >= SNIPPET_START_TIME) & (time_pad <= (SNIPPET_START_TIME + SNIPPET_DURATION));
plot(time_pad(peak_indices_within_snippet), signalC(peak_indices_within_snippet), 'ro')
hold off

time_snip = time_pad(snippet_start_index:snippet_end_index);
sigAsnip = signalA(snippet_start_index:snippet_end_index);
