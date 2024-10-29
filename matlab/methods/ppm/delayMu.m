function muListResult = delayMu(muListIn, m, truncate)
    arguments
        muListIn (:,:) double {mustBeInteger, mustBeNonnegative}
        m (1,1) double {mustBeInteger, mustBeNonnegative} = 1
        truncate logical = true
    end
    
    [countMuList, M] = size(muListIn);

    countZerosToPrepend = m;
    if truncate && countZerosToPrepend > M
        countZerosToPrepend = M;
    end
    
    if truncate
        countMuToKeep = max([M-m, 0]);
    else
        countMuToKeep = M;
    end

    zerosToPrepend = zeros(countMuList, countZerosToPrepend);
    muToKeep = muListIn(:, 1:countMuToKeep);
    
    muListResult = [zerosToPrepend, muToKeep];
end