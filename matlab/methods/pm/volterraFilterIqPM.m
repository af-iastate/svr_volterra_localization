function y = volterraFilterIqPM(x, h, M, K)
    N = length(x);
    M_iq = 2 * M;
    nSize = sum(arrayfun(@(k) nchoosek(M_iq+k-1, k), 1:K));
    xCompanion = zeros(nSize, 1);
    xSliceI = zeros(M, 1);
    xSliceQ = zeros(M, 1);
    y = zeros(N, 1);
    mtable = makeTablePM(M_iq, K);

    for ii=1:N
        xSliceI = [real(x(ii)); xSliceI(1:end-1)];
        xSliceQ = [imag(x(ii)); xSliceQ(1:end-1)];
        xCompanion = applyTableKernelPM(xCompanion, [xSliceI; xSliceQ], M_iq, mtable);
        y(ii) = xCompanion.' * h;
    end

end