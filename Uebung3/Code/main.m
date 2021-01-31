% Physikalische Geodaesie Uebung 3
% Nicholas Schneider & Ziqing Yu
% 27/01/2021

%% Initial
clc
close all
clearvars

%% Aufgabe 2 Teil a
% Implementierung und Visualisierung
phi = 0:pi/180:pi;
St = 1 ./ sin(phi/2) - 6 * sin(phi/2) + 1 - 5 * cos(phi) - 3 * cos(phi) .* log(sin(phi/2) + (sin(phi/2)).^2);
f = figure;
plot(phi/pi*180,St)
hold on
plot(phi/pi*180,zeros(length(0:180)));
xlabel('\psi [Grad]'); ylabel('St(\psi)');
legend('Stokes Funktion','Null')
saveas(f,'stokes','png');

% Nullstellen
syms phis
% Einheit: Grad
phi01 = double(vpasolve(1 / sind(phis/2) - 6 * sind(phis/2) + 1 - 5 * cosd(phis) - 3 * cosd(phis) * log(sind(phis/2) + (sind(phis/2))^2),phis,[30,50]));
phi02 = double(vpasolve(1 / sind(phis/2) - 6 * sind(phis/2) + 1 - 5 * cosd(phis) - 3 * cosd(phis) * log(sind(phis/2) + (sind(phis/2))^2),phis,[110,130]));

% Umrechnen in Grad, Minuten und Sekunden
[g1,m1,s1] = g2gms(phi01);
[g2,m2,s2] = g2gms(phi02);
% ganze Sekunde 
s1 = floor(s1);
s2 = floor(s2);

%% Aufgabe 2 Teil b

% Schwer Anomalie
ga = importdata('gravity_anomalies.txt');
ga_all = ga(:,[1,3,5,7,9]);

% Koordinaten
P1 = [48.40067893, 9.97228199]; % [phi,lambda] in Grad
P2 = [48.70311236, 9.65402314];
P3 = [48.80556353, 9.21339955];

N1 = GeoidHeight(P1,ga_all,phi01,phi02);
N2 = GeoidHeight(P2,ga_all,phi01,phi02);
N3 = GeoidHeight(P3,ga_all,phi01,phi02);


function[g,m,s]=g2gms(r)
% Umrechnen von dezimale Grad in Grad, Minuten und Sekunden
% Eingabe: r in Dezimal Grad 
% Ausgabe : g:Grad
%           m:Minute
%           s:Sekunde
if r>0
    g=floor(r);
    m=floor((r-g)*60);
    s=(r-g-m/60)*3600;
else
    r=-r;
    g=floor(r);
    m=floor((r-g)*60);
    s=(r-g-m/60)*3600;
    g=-g;
    m=-m;
    s=-s;
end
end

function[SA] = SphAbs(lambda1,lambda2,phi1,phi2)
% Eingabe: lambda1, lambda2, phi1, phi2: Länge und Breite von Beiden
% Punkten (in Grad)

% Ausgabe: sphärischer Abstand auf der Einheitskugel (in Grad)

SA = acosd(sind(phi1) .* sind(phi2) + cosd(phi1) .* cosd(phi2) .* cosd(lambda1 - lambda2));
end

function[StFun] = StFun(phi)
% phi in Grad
StFun = 1 ./ sind(phi/2) - 6 * sind(phi/2) + 1 - 5 * cosd(phi) - 3 * cosd(phi) .* log(sind(phi/2) + (sind(phi/2)).^2);
end

function [N] = GeoidHeight(P,ga_all,phi01,phi02)
R = 6371000;  % m
gamma = 9.81; % m/s^2 

Bmax = ga_all(:,1);
Bmin = ga_all(:,2);
Lmax = ga_all(:,3);
Lmin = ga_all(:,4);

PM = [(Bmax + Bmin)/2, (Lmax + Lmin)/2]; % Koordinaten der Mittelpunkten in Grad
dg = ga_all(:,5) * 1e-5; % Schwer Anomalie in m/s^2

Phi = SphAbs(P(2),PM(:,2),P(1),PM(:,1));  % in Grad
St_F = StFun(Phi);
Fla_Int = (Lmax - Lmin)/180*pi .* (sind(Bmax) - sind(Bmin)); 

id = (Phi~=phi01 | Phi~=phi02);

St_F = St_F(id);
dg = dg(id);
Fla_Int = Fla_Int(id);

dN = R / (4 * gamma * pi) .* St_F .* dg .* Fla_Int;
N = sum(dN);
end