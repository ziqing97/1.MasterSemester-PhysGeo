function[binom] = binomfun(a,b)
    binom = 1;
    for i = 1:b
        binom = binom * (a+1-i) / i;
    end
end