fun = @(b,c) tempo(b,c);
h1 = quad2d(fun,0,1,0,1);


function [a] = tempo(b,c)
    a = b ./ c;
end