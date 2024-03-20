function [allPeaks, allIndices, mad, filteredOutMad, medianY, filteredOutMedianY, firstThresholdY, secondThresholdY] = processChunks(arrValues, arrIndexes, highAmpFilt, transientsThresh, MIN_PK_WIDTH)

    allPeaks = {};
    allIndices = {};

    for i = 1:length(arrValues)
        chunk = arrValues(:,i);
        chunkIndexes = arrIndexes(:,i);

        chunk = chunk(~isnan(chunk));
        median = nanmedian(chunk);

        mad = nanmedian(abs(chunk-median));

        firstThreshold = median + (highAmpFilt*mad);

        greaterThanMad = find(chunk>firstThreshold);

        arr = 1:length(chunk);
        lowerThanMad = ~ismember(arr, greaterThanMad);
        filteredOut = chunk(find(lowerThanMad==1));

        filteredOutMedian = nanmedian(filteredOut);
        filteredOutMad = nanmedian(abs(filteredOut-nanmedian(filteredOut)));
        secondThreshold = filteredOutMedian+(transientsThresh*filteredOutMad);

        greaterThanThreshIndex = find(chunk>secondThreshold);
        greaterThanThreshValues = chunk(greaterThanThreshIndex);
        temp = zeros(1, length(chunk));
        temp(greaterThanThreshIndex) = greaterThanThreshValues;
        [~, peaks] = findpeaks(temp,'MinPeakDistance', MIN_PK_WIDTH);

        firstThresholdY = ones(1, length(chunk)) * firstThreshold;
        secondThresholdY = ones(1, length(chunk)) * secondThreshold;

        newPeaks = nan(1, length(chunk));
        newPeaks(peaks) = peaks + chunkIndexes(1);

        medianY = ones(1, length(chunk)) * median;
        filteredOutMedianY = ones(1, length(chunk)) * filteredOutMedian;

        allPeaks{end+1} = newPeaks;
        allIndices{end+1} = peaks;
    end
end