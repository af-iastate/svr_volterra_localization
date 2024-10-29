function y = volterraFilterIqDirect(x, h, M, K)
    N = length(x);
    M_iq = 2 * M;
    nSize = sum(arrayfun(@(k) nchoosek(M_iq+k-1, k), 1:K));
    xSlice = zeros(M, 1);
    y = zeros(N, 1);

    for ii=1:N
        xSlice = [x(ii); xSlice(1:end-1)];
        xCompanion = voltVecGen([real(xSlice); imag(xSlice)], 1:K);
        y(ii) = xCompanion.' * h;
    end

end