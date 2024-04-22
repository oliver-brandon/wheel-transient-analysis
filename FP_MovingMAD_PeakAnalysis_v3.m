clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = "3";
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
HAFT_THRESH_SNR = 'false';
includeHAFT = 'true';
MAD_MULTIPLIER = 3;
MIN_PK_WIDTH = 0.2 * SAMPLE_RATE;
%Snippet Args%
SNIPPET_DURATION = 300; % Duration of the snippet in seconds
SNIPPET_START_TIME = 1025; % Start time of the snippet in seconds
%Stream Stores%
DLS_ISOS = 'x405A'; % name of the 405A store
DLS_DA = 'x465A'; % name of the 465A store
NAC_ISOS = 'x405C'; % name of the 405C store
NAC_DA = 'x465C'; % name of the 465C store

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Version: %s\n",VERSION)
myDir = uigetdir('/Users/brandon/personal-drive/wheel-peak-test/tanks');
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

if strcmp(HAFT_THRESH_SNR, 'true')
    % Calculate signal to noise ratio for HAFT %
    %RMS%
    VnoiseA = rms(data.streams.(DLS_ISOS).data, 'omitnan');
    VsignalA = rms(data.streams.(DLS_DA).data, 'omitnan');
    VnoiseC = rms(data.streams.(NAC_ISOS).data, 'omitnan');
    VsignalC = rms(data.streams.(NAC_DA).data, 'omitnan');

    HAFT_THRESHA = VsignalA/VnoiseA;
    HAFT_THRESHC = VsignalC/VnoiseC;
else
    HAFT_THRESHA = 2;
    HAFT_THRESHC = 2;
end

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
signal_chunksA = [];
signal_chunksC = [];
signal_chunksAind = [];
signal_chunksCind = [];
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

% fills signal_chunksA and signal_chunksC with the actual signal chunks
for ii = 1:window_size:length(signalA)
    x1 = signalA(1,ii:ii+window_size-1);
    x2 = signalA_ind(1,ii:ii+window_size-1);
    signal_chunksA = [signal_chunksA; x1];
    signal_chunksAind = [signal_chunksAind; x2];
end
for ii = 1:window_size:length(signalC)
    x1 = signalC(1,ii:ii+window_size-1);
    x2 = signalC_ind(1,ii:ii+window_size-1);
    signal_chunksC = [signal_chunksC; x1];
    signal_chunksCind = [signal_chunksCind; x2];
end

% finds peaks within each chunk of signalA
for i = 1:height(signal_chunksA)
    
    [allPeaksA, allIndicesA, mad, filteredOutMad, medianY, filteredOutMedianY, ...
        firstThresholdY, secondThresholdY] = processChunks(signal_chunksA, ...
        signal_chunksAind, HAFT_THRESHA, MAD_MULTIPLIER, MIN_PK_WIDTH, includeHAFT);
end
% First, calculate the starting index of each chunk in the original signal
chunk_starts = (0:size(signal_chunksA, 2)-1) * size(signal_chunksA, 1);
% Adjust the peak indices to make them relative to the original signal
adjusted_peaksA = cellfun(@(peaks, start) peaks + start, allIndicesA, num2cell(chunk_starts), 'UniformOutput', false);
% Concatenate the cell arrays of adjusted peak indices and chunks to get a 1-D array of peaks and the original signal
adjusted_peaksA = cell2mat(adjusted_peaksA);
% Plots peaks for signalA
visualize_peaks(myDir, signalA, time_pad, adjusted_peaksA);

% finds peaks within each chunk of signalC
for i = 1:height(signal_chunksC)
    
    [allPeaksC, allIndicesC, mad, filteredOutMad, medianY, filteredOutMedianY, ...
        firstThresholdY, secondThresholdY] = processChunks(signal_chunksC, ...
        signal_chunksCind, HAFT_THRESHC, MAD_MULTIPLIER, MIN_PK_WIDTH, includeHAFT);
end
% First, calculate the starting index of each chunk in the original signal
chunk_starts = (0:size(signal_chunksC, 2)-1) * size(signal_chunksC, 1);
% Adjust the peak indices to make them relative to the original signal
adjusted_peaksC = cellfun(@(peaks, start) peaks + start, allIndicesC, num2cell(chunk_starts), 'UniformOutput', false);
% Concatenate the cell arrays of adjusted peak indices and chunks to get a 1-D array of peaks and the original signal
adjusted_peaksC = cell2mat(adjusted_peaksC);
% Plots peaks for signalC
visualize_peaks(myDir, signalC, time_pad, adjusted_peaksC);

% Count the number of peaks for signalA and signalC
numPeaksA = length(adjusted_peaksA)
numPeaksC = length(adjusted_peaksC)

% Calculate the frequency of peaks using SESSION_DURATION for time
frequencyA = numPeaksA / SESSION_DURATION
frequencyC = numPeaksC / SESSION_DURATION