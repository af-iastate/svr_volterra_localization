function [y, tDuration] = volterraFilterPM(x, h, M, K)
    N = length(x);
    nSize = sum(arrayfun(@(k) nchoosek(M+k-1, k), 1:K));
    xCompanion = zeros(nSize, 1);
    xSlice = zeros(M, 1);
    y = zeros(N, 1);
    mtable = makeTablePM(M, K);

    myTimer = tic;
    for ii=1:N
        xSlice = [x(ii); xSlice(1:end-1)];
        xCompanion = applyTableKernelPM(xCompanion, xSlice, M, mtable);
        y(ii) = xCompanion.' * h;
    end
    timeElapsed = toc(myTimer);

    if nargout > 1
        tDuration = timeElapsed;
    end
end