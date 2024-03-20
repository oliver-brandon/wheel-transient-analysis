function visualize_peaks(myDir, signal, time, peaksIndex)
    % Get the directory name and base name of the file
    [~, basename, ~] = fileparts(myDir);
    
    % Create a new figure
    figure;
    
    % Plot the z_score
    plot(time, signal, '-');
    xlim([0 3600])
    hold on;
    
    % Plot the peaks
    plot(time(peaksIndex), signal(peaksIndex), 'o');
    
    % Set the title of the plot
    title(basename);
    
    % Show the plot
    hold off;
end