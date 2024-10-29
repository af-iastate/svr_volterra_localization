function [xv] = voltVecGen(x,orders,fnc)
    % Andrew Bolstad, Oct. 8, 2019
    Nord = length(orders);
    xv = [];
    
    if nargin<3
        fnc = @recursiveTriangularForm;
    end
    
    for c = 1:Nord
        y = feval(fnc,x,orders(c));
        xv = [xv; y];
    end
end

