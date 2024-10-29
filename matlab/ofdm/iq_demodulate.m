function varargout = iq_demodulate( ...
        rx_in, Fs, fc, BW, options)
    arguments
        rx_in (:, 1) double
        Fs double {mustBePositive}
        fc double {mustBePositive}
        BW double {mustBePositive}
        options.phi double = 0
        options.filter_order double {mustBeInteger, mustBePositive} = 5
    end
    phi = options.phi;
    filter_order = options.filter_order;

    Ts = 1/Fs;
    t = Ts.*(0:length(rx_in)-1).';

    % (headroom_coeff*scs*nFFT/2)
    [B, A] = butter(filter_order, BW / (Fs/2), 'low');

    rx_i = 2*filtfilt(B, A, rx_in .*  cos(2*pi*fc.*t + phi));
    rx_q = 2*filtfilt(B, A, rx_in .* -sin(2*pi*fc.*t + phi));
    rx_iq = rx_i + 1j.*rx_q;

    if nargout >= 1
        varargout{1} = rx_iq;
    end
end