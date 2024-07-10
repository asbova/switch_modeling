function neurons = clusterStriatalNeurons(neurons, probabilityThreshold)

% Description
%
% Input:
%   neurons:    Structure containing all of the neurons for a dataset
%
% Output:

    % Constants
    options = statset('Display', 'final');
    cmap = lines(4);

    allWaveformData = NaN(length(neurons), 2);
    for iNeuron = 1 : length(neurons)
        allWaveformData(iNeuron, 1) = neurons(iNeuron).waveform.halfPeakWidth;
        allWaveformData(iNeuron, 2) = neurons(iNeuron).waveform.peakToTroughDuration;
    end

    selectUnits = allWaveformData(:,1) < 0.4 & allWaveformData(:,1) ~= 0; % Exclude outliers.
    selectWaveformData = allWaveformData(selectUnits, :);

    % Use a gaussian mixture model to separate neurons into clusters.
    gmModel = fitgmdist(selectWaveformData, 2, 'Options', options);
    neuronGroups = cluster(gmModel, selectWaveformData);
    neuronProbability = posterior(gmModel, selectWaveformData);

    % Find neurons that have probability of belonging to a cluster of > threshold (defined at top).
    probabilityIndex = NaN(size(neuronGroups));
    probabilityIndex(neuronProbability(:,1) > probabilityThreshold) = 1;
    probabilityIndex(neuronProbability(:,2) > probabilityThreshold) = 2;

    % Interneurons should have smaller peak to trough durations.
    if mean(selectWaveformData(probabilityIndex == 1, 2)) < mean(selectWaveformData(probabilityIndex == 2, 2))
        interneuronIndex = (probabilityIndex == 2);
        nonInterneuronIndex = (probabilityIndex == 1);
        probabilityIndex(interneuronIndex) = 1;
        probabilityIndex(nonInterneuronIndex) = 2;
    end

    % Identify MSNs vs. Interneurons vs. non-neurons
    msnCluster = (probabilityIndex == 1);
    interneuronCluster = (probabilityIndex == 2);   
    selectIndex = find(selectUnits);
    for iNeuron = 1 : length(selectIndex)
        if probabilityIndex(iNeuron) == 1
            neurons(selectIndex(iNeuron)).type = 'MSN';
        elseif probabilityIndex(iNeuron) == 2
            neurons(selectIndex(iNeuron)).type = 'INT';
        else
            neurons(selectIndex(iNeuron)).type = 'NA';
        end
    end

    % Output results
    fprintf('\nTotal MSNs: %d', sum(probabilityIndex == 1));
    fprintf('\nTotal Interneurons: %d', sum(probabilityIndex == 2));
    fprintf('\nTotal Excluded Neurons: %d\n', sum(isnan(probabilityIndex)));

    figure(1);
    % Plot all neurons: half peak width vs. peak to trough.
    subplot(1,2,1); cla;
    plot(allWaveformData(:,1), allWaveformData(:,2), '*');
    xlim([0 0.35])
    xlabel('Half Peak Width (ms)');
    ylabel('Peak To Trough (ms)');
    title('All Units')
    
    subplot(1,2,2); cla;
    hold on;
    scatter(selectWaveformData(msnCluster, 1), selectWaveformData(msnCluster, 2), 600, cmap(3,:), '.');
    scatter(selectWaveformData(interneuronCluster, 1), selectWaveformData(interneuronCluster, 2), 600, cmap(2,:), '.');
    xlim([0 0.35])
    xlabel('Half Peak Width (ms)');
    ylabel('Peak To Trough (ms)');
    text(0.025, 0.95, sprintf('MSNs: %d', sum(probabilityIndex == 1)), 'Color', cmap(3,:));
    text(0.025, 0.925, sprintf('INTs: %d', sum(probabilityIndex == 2)), 'Color', cmap(2,:));
    title('Clustered Units')

    saveas(gcf, fullfile('./results', sprintf('clustering_%0.2f.png', probabilityThreshold)));
    close all;