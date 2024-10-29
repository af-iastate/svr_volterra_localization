function [data, ndivs] = genPulsedToneSeries(freqs, duration, varargin)
    p = inputParser;
    p.addRequired('freqs');
    p.addRequired('duration');
    p.addOptional('backoff', []);
    p.addOptional('fs', 125e6);
    p.addOptional('npdFunction', @(x) x);
    p.addOptional('initialPadding', 12);
    p.addOptional('dutyCycle', 0.5);
    p.parse(freqs, duration, varargin{:});
    
    % hard-coded parameters
    memChunk = 1024;
    avePowerSetPoint = 0.28; % for 0 dB backoff
    
    % parameter init
    freqs = p.Results.freqs;
    duration = p.Results.duration;
    backoff = p.Results.backoff;
    fs = p.Results.fs;
    npdFunction = p.Results.npdFunction;
    initialPadding = p.Results.initialPadding;
    invDutyCycle = 1/(p.Results.dutyCycle);
    
    % freqs
    [nPulses, nTones] = size(freqs);
    
    % backoff
    if isempty(backoff)
        backoff = 6*ones(nPulses,1);
    elseif length(backoff)==1
        backoff = backoff*ones(nPulses,1);
    end
    nPulses = min(nPulses, length(backoff));
    
    % assume duration is in ms; add check later
    if length(duration)==1
        duration = [duration, duration/2];
    end
    
    % calculated params
    Ts = 1/fs;
    nSamples = floor(duration/Ts);
    p2pSamples = round(invDutyCycle*nSamples(1));
    if nargout>1
        ndivs = p2pSamples;
    end
    
    %Put tones on bins
    fres = fs/nSamples(2);
    toneBins = round(freqs/fres)+1;

    % average power is: nTones*A^2/2
    baseline = 10*log10(2*avePowerSetPoint/nTones);
    Adb = baseline - backoff;
    
    % make signal --------------------------------------------------------
    % loop over pulses
    data = zeros(p2pSamples*nPulses,1);
    for c = 1:nPulses
        % create tones on frequncy bins and apply backoff
        d = zeros(nSamples(2), 1);
        d(toneBins(c,:)) = 1;
        d(nSamples(2) + 2-toneBins(c,:)) = 1;
        d = 10^(Adb(c)/20) * (nSamples(2)/2) * ifft(d);

        % apply npd to d -----------------
        % transient
        x_init = [zeros(initialPadding,1); d];
        d_init = npdFunction(x_init);
        d_init(1:initialPadding) = [];

        % steady state
        x_ss = [d(end - initialPadding + 1:end); d];
        d = npdFunction(x_ss);
        d(1:initialPadding) = [];

        % ignore turn off transient for now
        % --------------------------------

        % repeat tones to fill pulse
        reps = floor(nSamples(1) / nSamples(2));
        rm = rem(nSamples(1), nSamples(2));
        d = [d_init; repmat(d, reps - 1,1); d(1:rm)];

        % fill in one pulse of data vector
        data((c-1) * p2pSamples + 1 : (c-1) * p2pSamples + nSamples(1)) = d;
    end

    % zero pad
%     zp1 = 5*memChunk;
%     data = [zeros(zp1,1); data; zeros(zp1,1)];
%     N = ceil(length(data)/memChunk)*memChunk;
%     zp2 = N-length(data);
%     data = [data; zeros(zp2,1)];
end