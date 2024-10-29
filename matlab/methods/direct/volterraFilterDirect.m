function [y, tDuration] = volterraFilterDirect(x, h, M, K)
    N = length(x);
    xSlice = zeros(M, 1);
    y = zeros(N, 1);
    myTimer = tic;
    for ii=1:N
        xSlice = [x(ii); xSlice(1:end-1)];
        xCompanion = voltVecGen(xSlice, 1:K);
        y(ii) = xCompanion.' * h;
    end
    timeElapsed = toc(myTimer);

    if nargout > 1
        tDuration = timeElapsed;
    end
end