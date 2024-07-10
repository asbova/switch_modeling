function neuronData = extractNeuronData(sessionData)

% Description
%
% Input:
%   sessionData:    Structure containing 
%
% Output:

    % Constants
    neuronData = struct;
    unitLabels = {'a', 'b', 'c', 'd', 'e', 'f'};
    upsamplingFactor = 6;
    sampleRate = 60000;
    timeUnitConversion = 1000; % Seconds to milliseconds.

    pl2Filename = sessionData.pl2FilePathway;
    pl2Info = PL2GetFileIndex(pl2Filename);

    % Find how many units are in each channel.
    nSpikeChannels = numel(pl2Info.SpikeChannels);
    nUnits = zeros(nSpikeChannels, 1);
    for iChannel = 1 : nSpikeChannels
        nUnits(iChannel) = sum(pl2Info.SpikeChannels{iChannel}.UnitCounts(2:end) > 0);
    end
    channelsWithUnits = find(nUnits > 0);

    unitIndex = 1;
    for iChannel = 1 : length(channelsWithUnits)
        currentChannel = channelsWithUnits(iChannel);
        for jUnit = 1 : nUnits(currentChannel)
            % Extract neuron waveforms and spikes.
            unitName = sprintf('%s%s', pl2Info.SpikeChannels{currentChannel}.Name, unitLabels{jUnit});
            waveformData = internalPL2Waves(pl2Filename, pl2Info.SpikeChannels{currentChannel}.Name, jUnit);
            waveformMean = mean(waveformData.Waves, 1);
            [~, spikeTimestamps] = plx_ts(pl2Filename, pl2Info.SpikeChannels{currentChannel}.Name, jUnit);
            
            neuronData(unitIndex).name = unitName;
            neuronData(unitIndex).filename = sessionData.pl2Name;
            neuronData(unitIndex).filePathway = pl2Filename;
            neuronData(unitIndex).spikeTimestamps = spikeTimestamps;
            neuronData(unitIndex).events = sessionData.ephysEvents;
            neuronData(unitIndex).waveform.mean = waveformMean;
            neuronData(unitIndex).channelNames = pl2Info.SpikeChannels{currentChannel}.Name;
            neuronData(unitIndex).unitNumber = jUnit;
            neuronData(unitIndex).nTimestamps = length(neuronData(unitIndex).spikeTimestamps);
            neuronData(unitIndex).averageFiringRate = neuronData(unitIndex).nTimestamps./max(neuronData(unitIndex).spikeTimestamps);

            unitIndex = unitIndex + 1;
        end
    end

    for iUnit = 1 : length(neuronData)
        upsampledWaveform = interp(neuronData(iUnit).waveform.mean, upsamplingFactor); % Upsample to 40*4 KHz.
        
        % Calculate the half peak width.
        [peakAmplitude, peakTime] = min(upsampledWaveform);
        [halfPeakAmplitude1, halfPeakTime1] = min(abs(upsampledWaveform(1:peakTime) - peakAmplitude/2));
        [halfPeakAmplitude2, halfPeakTime2] = min(abs(upsampledWaveform(peakTime:end) - peakAmplitude/2));
        halfPeakTime2 = halfPeakTime2 + peakTime - 1;
        halfPeakWidth = abs(halfPeakTime2 - halfPeakTime1)/sampleRate/upsamplingFactor*timeUnitConversion; % Milliseconds

        % Calculate the peak to trough.
        [troughAmplitude, troughTime] = max(upsampledWaveform(peakTime:end));
        troughTime = troughTime + peakTime - 1;
        peakToTroughDuration = (troughTime - peakTime)/sampleRate/upsamplingFactor*timeUnitConversion;
        peakToTroughRatio = abs(peakAmplitude/troughAmplitude);

        neuronData(iUnit).waveform.halfPeakWidth = halfPeakWidth;
        neuronData(iUnit).waveform.peakToTroughDuration = peakToTroughDuration;
        neuronData(iUnit).waveform.peakToTroughRation = peakToTroughRatio;
    end


end