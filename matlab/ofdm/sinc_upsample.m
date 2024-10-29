function xup = sinc_upsample(x, upfactor)
    nFFT = numel(x);
    nFFTUp = upfactor * nFFT;
    X = fft(x, nFFT);
    XUp = [
        X(1 : nFFT/2);
        zeros((upfactor-1) * nFFT, 1);
        X((nFFT/2 +1) : nFFT)
    ];
    
    % To time domain
    xup = upfactor * ifft(XUp, nFFTUp);
end

