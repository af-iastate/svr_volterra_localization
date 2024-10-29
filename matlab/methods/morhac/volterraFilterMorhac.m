function [y, tDuration] = volterraFilterMorhac(x, hMorhac, M, K)

    memories = repmat(M, [K, 1]);
    
    myTimer = tic;
    y = fastVMcell(x, hMorhac, memories);
    timeElapsed = toc(myTimer);
    
    y = sum(y, 1).';
    N = length(y);
    y = y(1:end - (N - length(x)));

    if nargout > 1
        tDuration = timeElapsed;
    end
end