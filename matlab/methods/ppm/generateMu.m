function muList = generateMu(M, p)
    arguments
        M (1,1) double {mustBeInteger, mustBePositive}
        p (1,:) double {mustBeInteger, mustBePositive}
    end
%GENERATEMU generate a matrix of multi-indices.
%   GENERATEMU(M, p) where M is the memory of the underlying indices
%   and p is a vector of the orders of the underlying indices.

    muList = cell([length(p), 1]);
    for ii=1:length(p)
        tauList = generateTauUpperTri(M, p(ii));
        [countTauList, ~] = size(tauList);
        
        % Convert taus -> mus
        muList{ii} = arrayfun(@(ii) ...
            arrayfun(@(x) sum(tauList(ii,:)==x), 0:M-1), ...
            1:countTauList, ...
            'UniformOutput', false ...
        );
        muList{ii} = cell2mat(muList{ii}.');
    end

    % Flatten
    muList = vertcat(muList{:});
    
end