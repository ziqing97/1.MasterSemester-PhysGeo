%% Physikalische Geodäsie Übung 1
% Aufgabe 3
% Nicholas Schneider & Ziqing Yu
% 30/11/2020
clc
clear all
close all

% 
% Drehwinkel
gamma = pi / 12;
% Laenge & Breite des Zentrums
alpha = (9 + 11/60) / 180 * pi;
beta = (48 + 46/60) /180 * pi;
beta = pi/2 - beta;

n = 6;
m = 4;

len_b = 180;
len_l = 360;

% P(6,4)
theta = linspace(0,pi,len_b);
P_all = Normalized_Lengendre(6,theta);
P_64 = P_all{7,5};
P_64 = repmat(P_64,len_l,1);

% Y(original)
lambda = linspace(-pi,pi,len_l);
lambda = lambda';
lambda = repmat(lambda,1,len_b);
lambda = lambda';
P_64 = P_64';
Y_64 = P_64 .* exp(1i * m *lambda) * (-1)^(m) / 2 / sqrt(2 * pi);
Y_64 = real(Y_64);
figure
hm = imagesc(Y_64);
colormap('jet')
title('Unrotiert')
ax = gca;
xticks(0:90:360);
xticklabels({'0','90','180','270','360'})
yticks(0:30:180);
yticklabels({'0','30','60','90','120','150','180'})
colorbar

% Y(gedreht)
Y_64d = 0;
for k = -n:n
    lambda2 = lambda - alpha;
    if (k<0)
        Pnk = P_all{7,(-k+1)};
        Pnk = Pnk';
        Ynk = Pnk .* exp(-1i * abs(k) *lambda2) / 2 / sqrt(2 * pi);
    elseif (k>0)
        Pnk = P_all{7,k+1};
        Pnk = Pnk';
        Ynk = Pnk .* exp(1i * abs(k) *lambda2) * (-1)^(abs(k)) / 2 / sqrt(2 * pi);
    else
        Pnk = P_all{7,1};
        Pnk = Pnk';
        Ynk = Pnk .* exp(1i * abs(k) *lambda2) / 2 / sqrt(pi);
    end
    Y_64d = Y_64d + Wigner_d(m,n,k,beta) * Ynk;
end

figure
Y_64d = Y_64d * exp(1i * m * gamma);
Y_64d = real(Y_64d);
hm2 = imagesc(Y_64d);
colormap('jet')
title('Rotation mit Wigner_d')
xticks(0:90:360);
xticklabels({'0','90','180','270','360'})
yticks(0:30:180);
yticklabels({'0','30','60','90','120','150','180'})
colorbar



