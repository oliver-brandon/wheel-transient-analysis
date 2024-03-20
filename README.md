# FP_MovingMAD_PeakAnalysis_v3.m

This MATLAB script is used for analyzing peaks in data streams. It was created by Brandon L. Oliver and adapted from Barker et al. (2017).

## Overview

The script starts by setting up various parameters such as the time threshold, downsampling rate, sample rate, session duration, window size, MAD multiplier, minimum peak width, snippet duration and start time, and the names of the data streams.

The script then prompts the user to select a directory containing the data to be analyzed. The data is loaded into the workspace using the TDTbin2mat function.

The time array for all streams is created and the initial part of the data, where the values may be erratic due to the LEDs being turned on, is removed.

The data streams and time array are then downsampled by a factor of N.

## Parameters

- `T`: Time threshold below which data will be discarded.
- `N`: Downsampling factor.
- `SAMPLE_RATE`: Sample rate of the data after downsampling.
- `SESSION_DURATION`: Duration of the session in seconds.
- `WINDOW_SIZE_SECONDS`: Size of the window for moving MAD calculation in seconds.
- `MAD_MULTIPLIER`: Multiplier for the MAD threshold.
- `MIN_PK_WIDTH`: Minimum width of a peak.
- `SNIPPET_DURATION`: Duration of the snippet in seconds.
- `SNIPPET_START_TIME`: Start time of the snippet in seconds.
- `DLS_ISOS`, `DLS_DA`, `NAC_ISOS`, `NAC_DA`: Names of the data streams.

## Usage

Run the script in MATLAB. You will be prompted to select a directory containing the data to be analyzed. The script will then process the data and perform peak analysis. 
