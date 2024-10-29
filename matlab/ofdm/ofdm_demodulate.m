function varargout = ofdm_demodulate( ...
        rx_iq, nFFT, options)
    arguments
        rx_iq (:, 1) double
        nFFT double {mustBePositive, mustBeInteger}
        options.upfactor double {mustBePositive, mustBeInteger} = 120
        options.nCycPre double {mustBePositive, mustBeInteger} = 32
    end
    
    upfactor = options.upfactor;
    nCycPre = options.nCycPre;
    
    nFFTUp = upfactor * nFFT;
    nCycPreUp = upfactor * nCycPre;
    
    % remove cyclic prefix
    rx_idxA = nCycPreUp+1;
    rx_sync = rx_iq(rx_idxA : end);

    % To freq domain
    rx_gridUp = fft(rx_sync(1:nFFTUp), nFFTUp) / upfactor;
    
    % downsample
    rx_grid = [
        rx_gridUp(1: nFFT/2);
        rx_gridUp((1+(upfactor-0.5)*nFFT) : end)
    ];

    % outputs
    if nargout >= 1
        varargout{1} = rx_grid;
    end
    if nargout > 1
        error('Too many output arguments. Expected at most 1. Found %d.', nargout)
    end
end



