function mapping = makeTablePPM(M, K)
    arguments
        M (1,1) double {mustBeInteger, mustBeNonnegative}
        K (1,1) double {mustBeInteger, mustBeNonnegative}
    end
    p = 1:K;
    muDest = generateMu(M, p);
    muSrc = delayMu(muDest, 1, false);
    % Increase the memory of the destination
    [N, ~] = size(muDest);
    muDest = [muDest, zeros(N, 1)];
    
    % mapping{1} = reassignMap, mapping{2} = permMap, mapping{3} = prodMap
    mapping = struct('reassign', [], 'permute', [], 'product', []);
    for ii=1:N
        muCurrent = muDest(ii,:);
        
        xNewPower = 0;
        previousTermIdx = 0;
        
        mostSignficant = muCurrent(1);
        if mostSignficant > 0
            muCurrent(1) = 0;
            xNewPower = mostSignficant;
        end
        % There is ABSOLUTELY a faster way to do this
        [wasFound, rowIdx] = ismember(muCurrent, muSrc, 'rows');
        altCheck = ~all(muCurrent == zeros(size(muCurrent)));
        if altCheck ~= wasFound
            error('\"wasFound\" (%d) differs from \"muCurrent == zeros\" (altCheck)', altCheck, wasFound)
        end
        if wasFound
            previousTermIdx = rowIdx;
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