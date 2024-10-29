function [freqsTrain, xTrain, freqsTest, xTest, splitIdx] = makeTonePulseDataSets(freqPool, durations, fs, numTones, numPulses, rngSeed)
    if nargin < 6
        rngSeed = [];
    end
    if ~isempty(rngSeed)
        rng(rngSeed);
    end
    
    numSetElements = numTones*numPulses;
    numDraws = 2 * numSetElements;
    
    freqsDrawn = freqPool(randperm(numel(freqPool), numDraws));
    idx = 1:numSetElements;
    freqsTrain = reshape(freqsDrawn(idx), numPulses, numTones);
    idx = idx + numSetElements;
    freqsTest = reshape(freqsDrawn(idx), numPulses, numTones);
    
    [xTrain, splitIdx] = getInputs(freqsTrain, durations, fs);
    xTest = getInputs(freqsTest, durations, fs);
end


function [x, splitIdx] = getInputs(freqsTwoTone, durationTwoToneAdditional, fs)
    [x, ndivs] = genPulsedToneSeries( ...
        freqsTwoTone, ...
        durationTwoToneAdditional, ...
        'dutyCycle', 0.5, ...
        'fs', fs ...
    );
    N = floor(length(x)/ndivs);
    splitIdx = (1:N) * ndivs;
end