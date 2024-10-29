function varargout = iq_modulate( ...
        tx_iq, Fs, fc, options)
    arguments
        tx_iq (:, 1) double
        Fs double {mustBePositive}
        fc double {mustBePositive}
        options.phi double = 0
    end
    phi = options.phi;
    Ts = 1/Fs;
    t = Ts.*(0:length(tx_iq)-1).';
    tx_out = real(tx_iq).*cos(2*pi*fc.*t + phi) - imag(tx_iq).*sin(2*pi*fc.*t + phi);

    if nargout >= 1
        varargout{1} = tx_out;
    end
end