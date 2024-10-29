function [xv] = voltVecGen(x,orders,fnc)
%voltVecGen Create polynomial combinations for Volterra filter
%   voltVecGen(x,orders) creates all the necessary polynomial combinations 
%   of the input x for implementing a Volterra filter with vector orders
%   specifying which order Volterra kernels to use.  The memory depth is
%   assumed to be length(x)-1.  It is assumed that x is of the form:
%           x = [ x(n) x(n-1) ... x(n-M) ].'
%
%   voltVecGen(x,orders,fnc) allows users to specify which Volterra model
%   to use via a function handle fnc.  By default, 
%   fnc = @recursiveTriangularForm()

% Andrew Bolstad, Oct. 8, 2019

%N = length(x);
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

