function varargout = ofdm_modulate( ...
        tx_grid, scs, options)
    arguments
        tx_grid (:, 1) double
        scs double {mustBePositive} = 15e3
        options.nCycPre double {mustBePositive, mustBeInteger} = 32
        options.upfactor double {mustBePositive, mustBeInteger} = 120
    end
    
    nFFT = numel(tx_grid);
    if mod(nFFT, 2) == 1
        error('Size of tx_grid must be even. Current size is: %d.', nFFT);
    end
    upfactor = options.upfactor;
    nCycPre = options.nCycPre;

    nFFTUp = upfactor * nFFT;
    nCycPreUp = upfactor * nCycPre;
    Fs = scs * nFFT/2;   % Base Sampling rate 
    FsUp = upfactor * Fs;

    % Oversample Grid
    tx_gridUp = [
        tx_grid(1 : nFFT/2);
        zeros((upfactor-1) * nFFT, 1);
        tx_grid((nFFT/2 +1) : nFFT)
    ];
    
    % To time domain
    tx_ofdm = upfactor * ifft(tx_gridUp, nFFTUp);
    tx_ofdm = tx_ofdm(:);

    % Add Cyclic Prefix
    tx_ofdmPrefix = tx_ofdm(nFFTUp-nCycPreUp+1 : nFFTUp);
    tx_ofdm = [tx_ofdmPrefix; tx_ofdm];

    % outputs
    if nargout >= 1
        varargout{1} = tx_ofdm;
    end
    if nargout >= 2
        varargout{2} = FsUp;
    end
    if nargout > 3
        error('Too many output arguments. Expected at most 2. Found %d.', nargout)
    end
end

