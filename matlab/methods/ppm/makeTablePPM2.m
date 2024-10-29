function mapping = makeTablePPM2(M, K)
    arguments
        M (1,1) double {mustBeInteger, mustBeNonnegative}
        K (1,1) double {mustBeInteger, mustBeNonnegative}
    end

    tauDest = generateTauUpperTriFull(M, K);
    N = size(tauDest, 1);

    % Precompute nchoosek
    nchoosek_table = zeros(M+K-1, K);
    for m=1:M+K-1
        for k=1:m
            nchoosek_table(m, k) = nchoosek(m, k);
        end
    end
    
    % mapping{1} = reassignMap, mapping{2} = permMap, mapping{3} = prodMap
    mapping = struct('reassign', [], 'permute', [], 'product', []);
    for ii=1:N
        
        xNewPower = 0;
        previousTermIdx = 0;
        
        tauCurrent = tauDest{ii};
        tauPrev = tauCurrent - 1;
        k = size(tauPrev, 2);
        xNewPower = sum(tauPrev < 0);
        if xNewPower < k
            previousTermIdx = tauToFlatIndex(M, tauPrev(xNewPower+1:end), 1:K, nchoosek_table) + 1;
        end

        if xNewPower > 0 && previousTermIdx > 0
            mapping.product(end+1,:) = [ ...
                ii, ...
                previousTermIdx, ...
                xNewPower ...
            ];

        elseif previousTermIdx > 0
            mapping.permute(end+1,:) = [ ...
                ii, ...
                previousTermIdx ...
            ];
        elseif xNewPower > 0
            mapping.reassign(end+1,:) = [ ...
                ii, ...
                xNewPower ...
            ];
        end
    end

end

function result = generateTauUpperTriFull(M, K)
        result = cell(K, 1);
        for k=1:K
            result{k} = generateTauUpperTri(M, k);
            result{k} = mat2cell(result{k}, ones(1, size(result{k}, 1)));
        end
        result = vertcat(result{:});
end

function result = tau2mu(tau, M)
    result = arrayfun(@(m) sum(tau==m), 0:M-1);
end

function result = mu2tau(mu, K)
    result = zeros(1, K);
    jj = 0;
    for ii=1:length(mu)
        amnt = mu(ii);
        for ll=1:amnt
            jj = jj + 1;
            result(jj) = ii-1;
        end
    end
    result = result(1:jj);
end