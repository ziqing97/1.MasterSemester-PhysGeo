function[binom] = binomfun(a,b)
% Diese Funktion berechnet die Binomialkoeffiyienten
    if (2*b)>a
        b = a - b;
    end
    binom = 1;
    for i = 1:b
        binom = binom * (a+1-i) / i;
    end
end