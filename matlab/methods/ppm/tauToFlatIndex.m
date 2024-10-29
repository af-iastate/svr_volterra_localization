function result = tauToFlatIndex(M, tau, p, nchoosek_table)
    
    function result = phi(M, tau, m)
        result = 0;
        result = result + sum(arrayfun(@(K) ...
            nchoosek_table(K+M-1, K), ...
            m ...
        ));
        result = result + tauToOffsetWithinOrder(M, tau, nchoosek_table); 
    end

    [nRows, K] = size(tau);
    m = p(p<K);
    
    if nRows == 1
        result = phi(M, tau, m);
        return
    end

    result = arrayfun(@(ii) phi(M, tau(ii, :), m), 1:nRows);
    result = result(:);
end