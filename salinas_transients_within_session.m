clear;
close all;
clc;
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dType = 2; % 1 = TDT, 2 = mat
t = 30; % first t seconds are discarded to remove LED on artifact
N = 100; % downsample signal N times
channel = 2; % 1 = A, 2 = C
ARTIFACT465 = 0;
saveArtifact = 0; % 0 = do not save, 1 = save, 2 = overwrite
figSave = 1; % saves to folder with mat files
prom = 3;
session_duration = 5400;
thresh_base1 = [t, 600]; % baseline for first 10 minutes (habituation)
thresh_base2 = [600, 1800]; % baseline for next 20 minutes (vehicle)
thresh_base3 = [1800, 5400]; % baseline for last 60 minutes (cocaine)
useOtherThresh = 1;
other_thresh1 = [t, 600];
other_thresh2 = [600, 1800];
other_thresh3 = [600, 1800];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dType == 1
    % Gets tank from UI pop-up window
    TANK_NAME = uigetdir('H:\\My Drive\\wheel-peak-test\\tanks', 'Select a tank to plot');
    if TANK_NAME == 0
        disp('Select a tank to start!')
        return
    end
    disp('Data Type: TDT Tank.')
    data = TDTbin2mat(TANK_NAME, 'T2', session_duration+t, 'TYPE', {'streams'});
elseif dType == 2
    [TANK_NAME, pathname] = uigetfile('E:\\Google Drive\\hot-wheels\\within-session-wheel-mats');
    
    if TANK_NAME == 0
        disp('Select a mat to start!')
        return
    end
    [~,name,~] = fileparts(TANK_NAME);
    TANK_NAME = strcat(pathname,TANK_NAME);
    load(TANK_NAME)
else
    disp('Select a data type')
end
if channel == 1
    ISOS = 'x405A';
    Grab = 'x465A';
elseif channel == 2
    ISOS = 'x405C';
    Grab = 'x465C';
else
    disp('Cannot resolve channel.')
    return
end

if strcmp(Grab,'x465A')
    ROI = 'DLS';
    thresh_mult = 0.5;
    
    if saveArtifact == 1 && ~isfield(data, "DLSartifact")
        data.DLSartifact = ARTIFACT465;
        disp('New artifact level saved')
    elseif saveArtifact == 2
        data.DLSartifact = ARTIFACT465;
        disp('Overwriting artifact level')
    elseif saveArtifact == 0 && isfield(data, "DLSartifact")
        ARTIFACT465 = data.DLSartifact;
        disp('Artifact level loaded')
    else
        disp('Artifact level not saved')
    end
elseif strcmp(Grab,'x465C') 
    ROI = 'NAc';
    thresh_mult = 0.3;
    
    if saveArtifact == 1 && ~isfield(data, "NACartifact")
        data.NACartifact = ARTIFACT465;
        disp('New artifact level saved')
    elseif saveArtifact == 2
        data.NACartifact = ARTIFACT465;
        disp('Overwriting artifact level')
    elseif saveArtifact == 0 && isfield(data, "NACartifact")
        ARTIFACT465 = data.NACartifact;
        disp('Artifact level loaded')
    else
        disp('Artifact level not saved')
    end
else
    disp('Cannot resolve ROI.')
    return
end

if saveArtifact == 1 || saveArtifact == 2
    disp('Saved.')
    save(TANK_NAME, 'data')
else
    disp('Artifact threshold not saved')
end

ISOS_raw = data.streams.(ISOS).data;
Grab_raw = data.streams.(Grab).data;

% Checks for unequal isosbestic and Grab sensor signal length correcting if
% necessary
if length(Grab_raw) < length(ISOS_raw)
    disp('Isosbestic signal array is longer than Grab signal array')
    ISOS_raw = ISOS_raw(1:length(Grab_raw));
    disp('Corrected.')
elseif length(Grab_raw) > length(ISOS_raw)
    disp('Isosbestic signal array is shorter than Grab signal array')
    Grab_raw = Grab_raw(1:length(ISOS_raw));
    disp('Corrected.')
else
    disp('Isosbestic and Grab signal arrays are of equal size')
    disp('No correction necessary.')
end

% time array
time = (1:length(Grab_raw))/data.streams.(Grab).fs;

% removes the first (t) seconds of signal
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
Grab_raw = Grab_raw(ind:end);
ISOS_raw = ISOS_raw(ind:end);

% Downsample streams and time array by N times%
ISOS_raw = downsample(ISOS_raw, N);
Grab_raw = downsample(Grab_raw, N);
time = downsample(time, N);



% dF/F
bls = polyfit(ISOS_raw, Grab_raw, 1);
Y_fit_all = bls(1) .* ISOS_raw + bls(2);
Y_dF_all = Grab_raw - Y_fit_all; %dF (units mV) is not dFF
Grab_dFF = 100*(Y_dF_all)./Y_fit_all;

% Photobleach correction
ISOS_raw = detrend(ISOS_raw);
Grab_dFF = detrend(Grab_dFF);

% noise reduction using moving median
Grab_filt = smoothdata(Grab_dFF,'movmedian',10);


% artifact removal
% Find indices where the signal is greater than the artifact threshold
art1_indices = find(any(Grab_filt > ARTIFACT465, 1));
% Find indices where the signal is less than the negative artifact threshold
art2_indices = find(any(Grab_filt < -ARTIFACT465, 1));
% Combine the indices of artifacts
artifact_indices = unique([art1_indices, art2_indices]);
% Create a logical array indicating good (non-artifact) signals
good_signals = true(1, size(Grab_filt, 2));
good_signals(artifact_indices) = false;
% Filter out the artifact signals
Grab_filt = Grab_filt(:, good_signals);
time_filt = time(:, good_signals);

f1 = figure(1);
% thresh_base1
subplot(5,1,1)
plot(time, Grab_dFF)
title('Without artifact removal')
ylabel('dFF')
xlim([t, session_duration])

subplot(5,1,2)
plot(time_filt, Grab_filt)
title('With artifact removal')
ylabel('dFF')
xlim([t, session_duration])

idx1 = find(thresh_base1(1) < time_filt & time_filt < thresh_base1(2));
if useOtherThresh == 0
    threshold1 = max(Grab_filt(:,idx1))*thresh_mult;
elseif useOtherThresh == 1
    alt_idx1 = find(other_thresh1(1) < time_filt & time_filt < other_thresh1(2));
    threshold1 = max(Grab_filt(:,other_thresh1))*thresh_mult;
end
[pks,locs] = findpeaks(Grab_filt(:,idx1), 'MinPeakHeight',threshold1, 'MinPeakProminence',prom);
numPeaks = length(locs);

peakFq = num2str((numPeaks/((thresh_base1(2) - thresh_base1(1) - t))*60), '%.4f');

subplot(5,1,3)
plot(time_filt(:, idx1), Grab_filt(:, idx1))
title(sprintf('ROI: %s, Baseline: %d, Prom: %d, ThreshMax: %.2f, Peaks: %d, Pk/m: %s', ROI, (thresh_base1(2) - thresh_base1(1)), prom, thresh_mult, numPeaks, peakFq))
ylabel('dFF')


xlim([thresh_base1(1), thresh_base1(2)]);

hold on
plot(time_filt(idx1(1) + locs - 1), pks, 'ro'); % Adjust locs by adding the starting index of the chunk
hold off


% thresh_base2
idx2 = find(thresh_base2(1) < time_filt & time_filt < thresh_base2(2));
if useOtherThresh == 0
    threshold2 = max(Grab_filt(:,idx2))*thresh_mult;
elseif useOtherThresh == 1
    alt_idx2 = find(other_thresh2(1) < time_filt & time_filt < other_thresh2(2));
    threshold2 = max(Grab_filt(:,other_thresh2))*thresh_mult;
end
[pks2,locs2] = findpeaks(Grab_filt(:,idx2), 'MinPeakHeight',threshold2, 'MinPeakProminence',prom);
numPeaks2 = length(locs2);
peakFq2 = num2str((numPeaks2/((thresh_base2(2) - thresh_base2(1)))*60), '%.4f');

subplot(5,1,4)
plot(time_filt(:,idx2), Grab_filt(:,idx2))
title(sprintf('ROI: %s, Baseline: %d, Prom: %d, ThreshMax: %.2f, Peaks: %d, Pk/m: %s', ROI, (thresh_base2(2) - thresh_base2(1)), prom, thresh_mult, numPeaks2, peakFq2))
ylabel('dFF')


xlim([thresh_base2(1), thresh_base2(2)]);

hold on
plot(time_filt(idx2(1) + locs2 - 1), pks2, 'ro'); % Adjust locs by adding the starting index of the chunk
hold off

% thresh_base3
idx3 = find(thresh_base3(1) < time_filt & time_filt < thresh_base3(2));
if useOtherThresh == 0
    threshold3 = max(Grab_filt(:,idx3))*thresh_mult;
elseif useOtherThresh == 1
    alt_idx3 = find(other_thresh3(1) < time_filt & time_filt < other_thresh3(2));
    threshold3 = max(Grab_filt(:,other_thresh3))*thresh_mult;
end
[pks3,locs3] = findpeaks(Grab_filt(:,idx3), 'MinPeakHeight',threshold3, 'MinPeakProminence',prom);
numPeaks3 = length(locs3);
peakFq3 = num2str((numPeaks3/(thresh_base3(2) - thresh_base3(1))*60), '%.4f');
subplot(5,1,5)
plot(time_filt(:,idx3), Grab_filt(:,idx3))
title(sprintf('ROI: %s, Baseline: %d, Prom: %d, ThreshMax: %.2f, Peaks: %d, Pk/m: %s', ROI, (thresh_base3(2) - thresh_base3(1)), prom, thresh_mult, numPeaks3, peakFq3))
ylabel('dFF')
xlabel('Time (s)')
xlim([thresh_base3(1), thresh_base3(2)]);
hold on
plot(time_filt(idx3(1) + locs3 - 1), pks3, 'ro'); % Adjust locs by adding the starting index of the chunk
hold off

% Display max dFF
fprintf('Max dFF: %.2f\n', max(Grab_dFF))

if figSave == 1
    disp('Saving figure.')
    % Save figure 1 to pathname
    file_name = char(strcat(pathname,name,{'_'},ROI,'.pdf'));
    orient(f1,'landscape');
    print(f1,file_name,'-dpdf','-vector','-bestfit','');
else
    disp('Figure not saved')
end

