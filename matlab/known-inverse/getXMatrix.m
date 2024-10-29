function X = getXMatrix(x, M, K, varargin)
% X = GETXMATRIX(x, M, p, dataRange)  Get volterra matrix from input x.
% Inputs:
%   x:         Input vector x.
%   M:         System memory M.
%   p:         Vector of term orders for which to compute the matrix.
%   dataRange: Range of x indices from which to construct X (default=M:end).
% Outputs:
%   X: Volterra matrix constructed from x
    
        % Parse arguments
        parser = inputParser;
        parser.addRequired('x');
        parser.addRequired('M');
        parser.addRequired('K');
        parser.addOptional('dataRange', []);
        parser.parse(x, M, K, varargin{:});
        dataRange = parser.Results.dataRange;

        % Calculate filter size
        nSize = sum(arrayfun(@(k) nchoosek(M+k-1, k), 1:K));

        % Process dataRange data
        N = length(x);
        if isempty(dataRange)
            dataRange = 1:N;
        elseif length(dataRange) == 1
            dataRange = dataRange:N;
        end

        
        xCompanion = zeros(1, nSize);
        X = zeros(length(dataRange), nSize);
        xSlice = zeros(M, 1);
        mtable = makeTablePM(M, K);
    
        for ii=1:N
            xSlice = [x(ii); xSlice(1:end-1)];
            xCompanion = applyTableKernelPM(xCompanion, xSlice, M, mtable);
            X(ii, :) = xCompanion(:).';
        end
        
end