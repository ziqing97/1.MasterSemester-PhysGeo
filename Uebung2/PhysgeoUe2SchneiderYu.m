clc
close all
clearvars
%%
% path zur shbundle, EGM96.mat war vorher erzeugt
addpath(genpath('E:\Studium\M1-Physgeo\PhysGeoUebung\shbundle-master\shbundle-master'))
load('PhysgeoUe2SchneiderYuEGM96.mat')

% Konstanten
GM = 3.986004415e14;
a = 6378136.701;
b = 6356751.661;
omega = 7.292115e-5;

r = 6371000; 
R = 6371000; 

l1 = 8;
l2 = 36;

% Breite und Co-breite
theta = linspace(0,pi,180);
phi = pi/2 - theta;

 
% L = 2,6,4,8
U_mittelwert = 0;
for i = 1:2
    U_mittelwert = U_mittelwert + (a/r)^(2*i) * J2l(i) .* Lengendpoly(2*i, theta);
end
U2 = GM / r * (1 - U_mittelwert);

U_mittelwert = 0;
for i = 1:4
    U_mittelwert = U_mittelwert + (a/r)^(2*i) * J2l(i) .* Lengendpoly(2*i, theta);
end
U4 = GM / r * (1 - U_mittelwert);

U_mittelwert = 0;
for i = 1:6
    U_mittelwert = U_mittelwert + (a/r)^(2*i) * J2l(i) .* Lengendpoly(2*i, theta);
end
U6 = GM / r * (1 - U_mittelwert);

U_mittelwert = 0;
for i = 1:8
    U_mittelwert = U_mittelwert + (a/r)^(2*i) * J2l(i) .* Lengendpoly(2*i, theta);
end
U8 = GM / r * (1 - U_mittelwert);

figure
hold on 
plot(U2)
plot(U4)
plot(U6)
plot(U8)

% Darstellung von Somigliana Pizzetti Normalpotential zum L=8
U_mittelwert = 0;
for i = 1:l1
    U_mittelwert = U_mittelwert + (a/r)^(2*i) * J2l(i) .* Lengendpoly(2*i, theta);
end
U8n = GM / r * (1 - U_mittelwert) + omega^2 / 2 * r^2 .* (cos(phi)).^2;
U8 = GM / r * (1 - U_mittelwert);
U8map = repmat(U8',1,360);

figure
hm = imagesc(U8map);
colormap('jet')
title('Somigliana-Pizzetti')
xticks(0:90:360);
xticklabels({'0','90','180','270','360'})
yticks(0:30:180);
yticklabels({'0','30','60','90','120','150','180'})
colorbar

% 'Volle' Schwerepotential mit EGM96
nm = EGM96(:,1:2)';
inC = EGM96(:,3)';
inS = EGM96(:,4)';
cs = vec2cs(nm, inC, inS);

a = 6378136.701;
b = 6356751.661;
omega = 7.292115e-5;
GM = 3.986004415e14;

f = gshs_grid(cs, linspace(0,2*pi,360), linspace(-pi/2,pi/2,180),r,'GM',GM, 'quant','potential','sub_wgs84',false);
figure
imagesc(f)
colormap('jet')
title('EGM96')
xticks(0:90:360);
xticklabels({'0','90','180','270','360'})
yticks(0:30:180);
yticklabels({'0','30','60','90','120','150','180'})
colorbar

% Differenz berechnen und darstellen
diff = f - U8map;
figure
imagesc(diff)
colormap('jet')
title('Diff')
xticks(0:90:360);
xticklabels({'0','90','180','270','360'})
yticks(0:30:180);
yticklabels({'0','30','60','90','120','150','180'})
colorbar

% Entlang Greenwich
Greenwich_voll = f(:,1);
Greenwich_sopi = U8map(:,1);
Greenwich_diff = Greenwich_voll - Greenwich_sopi;
figure
hold on 
plot(Greenwich_voll)
plot(Greenwich_sopi)
title('Potential entlang Greenwich')
legend('EGM96','Somigliana-Pizzetti')

figure
plot(Greenwich_diff)
title('Differenz entlang Greenwich')

function[J2l] = J2l(l)
% diese Funktion berechnet die J2l mit verschiedene l
GM = 3.986004415e14;
a = 6378136.701;
b = 6356751.661;
omega = 7.292115e-5; 

e1 = sqrt(a^2 - b^2)/a;
e2 = sqrt(a^2 - b^2)/b;
mhat = omega^2 * a^2 * b / GM;
q0 = 0.5 * ((3/e2^2 + 1) * atan(e2) - 3/e2);

J2 = e1^2/3 * (1 - 2 * e2 * mhat / 15 / q0);
if l == 1
    J2l = J2;
else
    J2l = (-1)^(l+1) * 3 * e1^(2*l) * (1 - l + 5 * l * J2 / e1^2) / (2 * l + 1) / (2 * l + 3);
end
end

function[Pl] = Lengendpoly(l,theta)
% diese Funktion berechnet die Lengendre Polynomen, die Eingabe sind l_max
% und Co-Breite. 
t = cos(theta);
P = cell(l,1);
P(:,:) = {NaN};
P(1) = {t};
if l > 1
    P(2) = {3/2 * t.^2 - 1/2}; 
end
if l == 1
    Pl = P{1};
elseif l == 2
    Pl = P{2};
else
    for i = 3:l
        P(i) = {(2 * i - 1) / l .* t .* P{i-1} - (l - 1) / l .* P{i-2}};
    end
    Pl = P{l};
end
end