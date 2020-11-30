function[d] = Wigner_d(m,n,k,t,beta)
% Nicholas Schneider & Ziqing Yu
% 30/11/2020
% Beschreibung der Funktion: 
% input:

% output:


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

% to be continued

end