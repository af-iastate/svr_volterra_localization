function [y, tDuration] = volterraFilterPPM(x, h, M, K)
    global volterraFilterPPM_makeTablePPM_memo;
    if isempty(volterraFilterPPM_makeTablePPM_memo)
        volterraFilterPPM_makeTablePPM_memo = memoize(@makeTablePPM2);
    end

    N = length(x);
    y = zeros(N, 1);

    mtable = volterraFilterPPM_makeTablePPM_memo(M, K);
    mreassign = mtable.reassign;
    mpermute = mtable.permute;
    mproduct = mtable.product;

    nSize = size(h, 1);
    xCompanion = zeros(nSize, 1);
    
    myTimer2 = tic;
    for ii=1:N
        xCompanion = applyTableKernelPPM(xCompanion, x(ii), K, ...
            mreassign, mpermute, mproduct);
        y(ii) = xCompanion.' * h;
    end
    timeElapsed = toc(myTimer2);

    if nargout > 1
        tDuration = timeElapsed;
    end
end