function h = GenBandstopIr(N, Fs, f0, f1, tp)
    if nargin<5
        tp = 'stop';
    end
    [B, A] = butter(1, [f0, f1]*2*pi, tp, 's');
    H = tf(B, A);
    G = c2d(H, 1/Fs, 'tustin');

    t = (0:N-1)/Fs;
    [h, ~] = impulse(G, t);
    h = h/Fs;
end