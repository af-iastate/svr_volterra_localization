function u = ThirdOrderActivationFunc(t, a1, a3)
    if a1 ~= 0 && a3 ~= 0
        tmpf = nthroot(-9*a3^2*t + sqrt(12*a1^3*a3^3 + 81*a3^4*t.^2), 3);
        u = (2/3)^(1/3)*a1./tmpf - tmpf./(2^(1/3)*3^(2/3)*a3);
    elseif a3 ~= 0
        u = nthroot(t./a3, 3);
    elseif a1 ~= 0
        u = (1/a1).*t;
    else
        u = 0.*t;
    end
end