function[delta_ij] = orthotest(i,j)
    if i == j
        delta_ij = 1;
    else
        delta_ij = 0;
    end
end