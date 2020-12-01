function[d] = Wigner_d(m,n,k,beta)
% Nicholas Schneider & Ziqing Yu
% 30/11/2020
% Mit dieser Funktion kann man die Wigner-d-Funktion berechnen. 
% input: Grad n,
%        Ordnung m und k
%        Winkel beta

% output: Wigner-d-Funktion


% pruefen, ob m,n,k ganzzahlig ist
if (rem(m,1)~=0 || rem(n,1)~=0 || rem(k,1)~=0)
    error('n,m,k müssen ganzzahlig sein.')
end
% pruefen, ob -n<=k und m<=n
if (-n > k || m> n)
    error('Bedingungen -n<=k oder m<=n nicht erfüllt.')
end
% SQRT berechnen
SQRT = sqrt((factorial(n+m) * factorial(n-m)) / (factorial(n+k) * factorial(n-k)));

%
t1 = max(0,k-m);
t2 = min(n-m,n+k);

Sigma = 0;
for t = t1:t2
    a = m - k + 2 * t;
    Sigma = Sigma + binomfun(n+k,t) * binomfun(n-k,n-m-t) * (-1)^t * cos(beta/2)^(2*n-a) * sin(beta/2)^a;
end
d = SQRT * Sigma;

end