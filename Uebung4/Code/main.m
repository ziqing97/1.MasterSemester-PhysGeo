%% PG Uebung 4
% Nicholas Schneider & Ziqing Yu
clc
close all
% clear all
% 
datafull = importdata('moon_Jan2000.txt');
data = datafull.data;
clear datafull

l = 2;
% Konstant
mM = 7.35e22; %kg
G = 6.672e-11; %m^3/kg*s^2
Re = 6378136.3; %m

r = data(:,2)*1000; % m
phi = data(:,4)/180*pi;  % radiant
lambda = data(:,3)/180*pi;  % radiant
theta = pi/2 - phi;

% a
Plm = Normalized_Lengendre(2,theta);
Plm2 = [Plm{3,3},Plm{3,2},Plm{3,1},Plm{3,2},Plm{3,3}];
faktor = [sin(2*lambda),sin(lambda),cos(0*lambda),cos(1*lambda),cos(2*lambda)];
Ylm2 = faktor.*Plm2;

vtidlm = G*mM./r /(2*l+1) .* (Re./r).^l .* Ylm2; 
vtidlm_15 = vtidlm(15,:);

% b
lambda_p = 8.33/180*pi; % rad
phi_p = 48.18/180*pi; % rad
theta_p = pi/2-phi_p; % rad
r_p = 6366837; % m

Plm_p = Normalized_Lengendre(2,theta_p);
Plm_p2 = [Plm_p{3,3},Plm_p{3,2},Plm_p{3,1},Plm_p{3,2},Plm_p{3,3}];
faktor_p = [sin(2*lambda_p),sin(lambda_p),cos(0*lambda_p),cos(1*lambda_p),cos(2*lambda_p)];
Ylm_p2 = Plm_p2.*faktor_p;
Ylm_p2 = repmat(Ylm_p2,20,1);

vtid = sum(vtidlm.*Ylm_p2.*(r_p/Re).^l,2);

vtid15 = vtid(1:5);

figure
plot(vtid)
xlabel('Tag')
ylabel('W')

% c
syms rs lambdas phis thetas
Plms(1) = sqrt(5)/4*(1-3*cos(2*phis));
Plms(2) = sqrt(15)/2*sin(2*phis); 
Plms(3) = sqrt(15)/4*(1+cos(2*phis));

clear Ylms
Ylms(1) = Plms(3)*sin(2*lambdas); % m=-2
Ylms(2) = Plms(2)*sin(lambdas);   % m=-1 
Ylms(3) = Plms(1);             % m=0
Ylms(4) = Plms(2)*cos(lambdas);   % m=1 
Ylms(5) = Plms(3)*cos(2*lambdas); % m=2


Ylms = repmat(Ylms,20,1);

vtids = sum(vtidlm .* Ylms*(rs/Re).^l,2);
vdrs = diff(vtids, rs);
vdlambdas = diff(vtids, lambdas);
vdphis = diff(vtids, phis);

gtids = [vdrs  vdlambdas/rs   vdphis/rs/cos(phi_p)];
gtid = eval(subs(gtids, [rs lambdas phis], [Re lambda_p phi_p]));


%% Aufgabe 2
t = [0:19]/36525;
t = t';

arg1 = 4.226311372109 + 221721.966199140683*t + 0.000030885540*t.^2 - 0.000000031964*t.^3 + 0.000000000154*t.^4;
arg1 = repmat(arg1,1,201);
% mean lunar longitude (s)  
arg2 = 3.810344278153 + 8399.709110962742*t - 0.000025593314*t.^2 + 0.000000032313*t.^3 - 0.000000000268*t.^4;
arg2 = repmat(arg2,1,201);
% mean solar longitude (h)
arg3 = 4.895062996673 + 628.331965369029*t + 0.000005292226*t.^2 + 0.000000000349*t.^3 - 0.000000000114*t.^4;
arg3 = repmat(arg3,1,201);
% mean longitude of lunar perigee (p)
arg4 = 1.454788534659 + 71.017685243656*t -0.000180148037*t.^2 -0.000000218021*t.^3 + 0.000000000919*t.^4;
arg4 = repmat(arg4,1,201);
% negative mean longitude of lunar ascending node (N_Strich)    
arg5 = 4.100746110564 + 33.757045953631*t -0.000036226248*t.^2 -0.000000037340*t.^3 + 0.000000000288*t.^4;
arg5 = repmat(arg5,1,201);
% mean longitude of solar perigee (ps)
arg6 = 4.938188176939 + 0.030010197632*t + 0.000007974215*t.^2 -0.000000000310*t.^3 -0.000000000058*t.^4;
arg6 = repmat(arg6,1,201);
% mean longitude of Mercury (LMer)
arg7 = 4.402608842461 + 2608.814705769090*t + 0.000005297047*t.^2 + 0.000000000316*t.^3 -0.000000000114*t.^4;
arg7 = repmat(arg7,1,201);
% mean longitude of Venus (LVen)
arg8 = 3.176146696956 + 1021.352941711328*t + 0.000005412955*t.^2 + 0.000000000260*t.^3 -0.000000000114*t.^4;
arg8 = repmat(arg8,1,201);
% mean longitude of Mars (LMar)
arg9 = 6.203476112911 + 334.085626078729*t + 0.000005419574*t.^2 + 0.000000000273*t.^3 - 0.000000000114*t.^4;
arg9 = repmat(arg9,1,201);
% mean longitude of Jupiter (LJup)
arg10 = 0.599547105074 + 52.993480508539*t + 0.000003897272*t.^2 + 0.000000000646*t.^3 -0.000000000091*t.^4;
arg10 = repmat(arg10,1,201);
% mean longitude of Saturn (LSat)
arg11 = 0.874016284019 + 21.354296582042*t + 0.000009059625*t.^2 - 0.000000000521*t.^3 -0.000000000170*t.^4;
arg11 = repmat(arg11,1,201);
fid = fopen('HW95_Katalog_Mond_threshold_Gr2.txt');
for i = 1:201
   	 number_thr(i) = fscanf(fid,'%i',1);
     planet1_thr(i)=fscanf(fid,'%1s',1);  
     planet2_thr(i)=fscanf(fid,'%1s',1);  
     degree_thr(i)=fscanf(fid,'%i',1);
     order_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k2_thr(:,i)=fscanf(fid,'%i',1) * ones(20,1);
     k3_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k4_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k5_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k6_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k7_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k8_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k9_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k10_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     k11_thr(:,i)=fscanf(fid,'%i',1)* ones(20,1);
     frequency_thr(i)=fscanf(fid,'%f',1);
     C0_thr(i)=fscanf(fid,'%f',1);      
   	 S0_thr(i)=fscanf(fid,'%f',1);  
   	 C1_thr(i)=fscanf(fid,'%f',1);
   	 S1_thr(i) =fscanf(fid,'%f',1); 
end       
fclose(fid);
for i = 1:20
    Clm(i,:) = C0_thr + t(i) * C1_thr;
    Slm(i,:) = S0_thr + t(i) * S1_thr;
end

syms phis lambdas rs
Plms(1) = sqrt(5)/4*(1-3*cos(2*phis));
Plms(2) = sqrt(15)/2*sin(2*phis); 
Plms(3) = sqrt(15)/4*(1+cos(2*phis));


alpha = order_thr*lambdas + order_thr.*arg1 + k2_thr.*arg2 + k3_thr.*arg3 +  k4_thr.*arg4...
    + k5_thr.*arg5 +  k6_thr.*arg6 +  k7_thr.*arg7 +  k8_thr.*arg8 +  k9_thr.*arg9 +  k10_thr.*arg10...
    + k11_thr.*arg11;

id0 = find(order_thr(i,:)==0);
id1 = find(order_thr(i,:)==1);
id2 = find(order_thr(i,:)==2);

for i = 1:20
    temp0 = sum(Clm(i,id0).*cos(alpha(i,id0)) + Slm(i,id0).*sin(alpha(i,id0)));
    temp1 = sum(Clm(i,id1).*cos(alpha(i,id1)) + Slm(i,id1).*sin(alpha(i,id1)));
    temp2 = sum(Clm(i,id2).*cos(alpha(i,id2)) + Slm(i,id2).*sin(alpha(i,id2)));
    VHW(i) = temp0 * Plms(1) * (rs/Re)^2 + temp1 * Plms(2) * (rs/Re)^2 + temp2 * Plms(3) * (rs/Re)^2;
    vdrsHW = diff(VHW(i), rs);
    vdlambdasHW = diff(VHW(i), lambdas);
    vdphisHW = diff(VHW(i), phis);
    gtidHWs(i,:) = [vdrsHW  vdlambdasHW/rs  vdphisHW/rs/cos(phis)];
end

gtidHW = eval(subs(gtidHWs, [rs lambdas phis], [Re lambda_p phi_p]));
