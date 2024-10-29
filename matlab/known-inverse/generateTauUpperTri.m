function result = generateTauUpperTri(M, K)
%GENERATETAUUPPERTRI(M,K) 
%Generates tau indices in upper-triangular order.
%   M: amount of memory (scalar)
%   K: order (scalar)
    
    function result = phi(offset, M, p)
        if length(p) == 1
            result = (offset:M-1).';
            return
        end
        
        resultCells = cell(M-offset, 1);
        for ii = offset:M-1
            innerResult = phi(ii, M, p(2:end));
            
            [N, ~] = size(innerResult);
            resultTmp = [ii*ones(N, 1), innerResult];
            resultCells{ii-offset+1} = resultTmp; 
        end
        
        result = vertcat(resultCells{:});
    end
    
    result = phi(0, M, 1:K);
end