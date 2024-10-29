function [y] = recursiveTriangularForm(x,p,M,xp,mp)
    % Andrew Bolstad, Oct. 15, 2019
    % initialize
    y = [];
    if nargin<3
        xp = 1;
        mp = 1;
        M = length(x);
        
        % assume user enters, e.g., p = 3 for third order
        p = 1:max(p);
    end
    
    if length(p)==1
        % innermost loop (really mp:M-1 with 1-index)
        y = xp*x(mp:M);
        y = y(:);
    else
        for m = mp:M
            % multiply by relevant x
            xcum = xp*x(m);
            
            % recurse
            ytmp = recursiveTriangularForm(x,p(2:end),M,xcum,m);
            
            y = [y; ytmp];
        end
    end
end

