function [hResult, KOut] = ApplyPolymap(varargin)
    [hIn, MIn, KIn, f, KOut] = parseLongForm(varargin{:});
    NPoly = length(f);
    hStructs = cell([NPoly, 1]);
    parfor q=1:NPoly
        [~, hStruct] = qthPower(hIn, MIn, KIn, q);
        hStructs{q} = hStruct;
    end
    
    KTotal = hStructs{NPoly}.K;
    hResult = [];
    idx = 1;
    
    if isempty(KOut)
        KOut = KTotal;
    end
    
    for ii=1:KOut
        tauArr = generateTauUpperTri(MIn, ii);
        [N, ~] = size(tauArr);
        multiIndices = arrayfun(@(ii) ...
            arrayfun(@(x) sum(tauArr(ii,:)==x), 0:MIn-1), ...
            1:N, ...
            'UniformOutput', false ...
        );
        multiIndices = cell2mat(multiIndices.');
    
        maxOrder = NPoly;
        minOrder = ceil(ii/KIn);
        for jj=1:N
            multiIndex = multiIndices(jj, :);
            res = 0;
            for kk=minOrder:maxOrder
                a = f(kk);
                res = res + a*hStructs{kk}.get(multiIndex);
            end
            hResult(idx) = res;
            idx = idx + 1;
        end
    end
    hResult = hResult(:);
end

function [h, M, K, f, KOut] = parseLongForm(varargin)
    parser = inputParser();
    parser.addRequired('h', @isnumeric);
    parser.addRequired('M', @isnumeric)
    parser.addRequired('K', @isnumeric)
    parser.addRequired('f', @isnumeric)
    parser.addOptional('KOut', [], @isnumeric)
    parser.parse(varargin{:});
    h = parser.Results.h;
    M = parser.Results.M;
    K = parser.Results.K;
    f = parser.Results.f;
    KOut = parser.Results.KOut;
end