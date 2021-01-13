clc
clearvars
close all

addpath(genpath('E:\Studium\M1-Physgeo\PhysGeoUebung\shbundle-master\shbundle-master'))
load('EGM96.mat')
nm = EGM96(:,1:2)';
inC = EGM96(:,3)';
inS = EGM96(:,4)';
cs = vec2cs(nm, inC, inS);

a = 6378136.701;
b = 6356751.661;
omega = 7.292115e-5;
GM = 3.986004415e14;

f = gshs_grid(cs, linspace(0,2*pi,360), linspace(-pi/2,pi/2,180),a,'GM',GM, 'quant','potential','sub_wgs84',false);
imagesc(f)
colormap('jet')
title('Somigliana-Pizzetti')
xticks(0:90:360);
xticklabels({'0','90','180','270','360'})
yticks(0:30:180);
yticklabels({'0','30','60','90','120','150','180'})
colorbar