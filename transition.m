function y = transition(x,spec)

arguments
    x (1,1) double {mustBeNonnegative}
    spec (1,:) char
end

if x >= 1
    y = 0;
else
    switch spec
        case 'indicator'
            y = 1;
        case 'mollifier'
            y = exp(1-1/(1-x^2));
        case 'cubic'
            y = 2*x^3 - 3*x^2 + 1;
        case 'linear'
            y = 1-x;
        case 'expReciprocal'
            y = exp(-1/(1-x)) / (exp(-1/x)+exp(-1/(1-x)));
    end
end