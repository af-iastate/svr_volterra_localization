function result = tauToOffsetWithinOrder(M, tau, nchoosek_table)
    
    function result = phi(M, tau, carryTerm)
        K = length(tau);
        if K == 1
            result = tau(1) - carryTerm;
            return;
        end

        subsetSize = @(M, K, tau0) nchoosek_table(M-tau0+K-1, K);
        
        result = subsetSize(M, K, carryTerm) - subsetSize(M, K, tau(1));
        result = result + phi(M, tau(2:end), tau(1));
    end
    
    [nRows, ~] = size(tau);
    if nRows == 1
        result = phi(M, tau, 0);
        return
    end
    
    result = arrayfun(@(ii) phi(M, tau(ii, :), 0), 1:nRows);
    result = result(:);
end