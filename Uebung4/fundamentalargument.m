function arg = fundamentalargument(t)
%  arg = fundamentalargument(t)
% 
% Berechnung der astronomischen Fundamentalargumente nach den Formeln und Bezeichnungen 
% ...
%
% EIN: 
% t :       die Zeit seit 1. Januar 2000, 12h UT in Julianischen Jahrhunderten
%           (1 Julianisches Jahrhundert = 36525 Tage
% 
% AUS:
% arg:     Cell-array, der 11 astronomischen Fundamentalargumente, jeweils in 
%          [rad]
%
% author:                   Hassan? /ANTONI
% project:                  Physkalische Geodäsie / Gezeitenübung
% created at:                   ??
% last modification:        01.02.2021 -> update to a function file

% mean local lunar time (tau) 
arg{1} = 4.226311372109 + 221721.966199140683*t + 0.000030885540*t.^2 - 0.000000031964*t.^3 + 0.000000000154*t.^4;

% mean lunar longitude (s)  
arg{2} = 3.810344278153 + 8399.709110962742*t - 0.000025593314*t.^2 + 0.000000032313*t.^3 - 0.000000000268*t.^4;

% mean solar longitude (h)
arg{3} = 4.895062996673 + 628.331965369029*t + 0.000005292226*t.^2 + 0.000000000349*t.^3 - 0.000000000114*t.^4;

% mean longitude of lunar perigee (p)
arg{4} = 1.454788534659 + 71.017685243656*t -0.000180148037*t.^2 -0.000000218021*t.^3 + 0.000000000919*t.^4;

% negative mean longitude of lunar ascending node (N_Strich)    
arg{5} = 4.100746110564 + 33.757045953631*t -0.000036226248*t.^2 -0.000000037340*t.^3 + 0.000000000288*t.^4;

% mean longitude of solar perigee (ps)
arg{6} = 4.938188176939 + 0.030010197632*t + 0.000007974215*t.^2 -0.000000000310*t.^3 -0.000000000058*t.^4;

% mean longitude of Mercury (LMer)
arg{7} = 4.402608842461 + 2608.814705769090*t + 0.000005297047*t.^2 + 0.000000000316*t.^3 -0.000000000114*t.^4;

% mean longitude of Venus (LVen)
arg{8} = 3.176146696956 + 1021.352941711328*t + 0.000005412955*t.^2 + 0.000000000260*t.^3 -0.000000000114*t.^4;

% mean longitude of Mars (LMar)
arg{9} = 6.203476112911 + 334.085626078729*t + 0.000005419574*t.^2 + 0.000000000273*t.^3 - 0.000000000114*t.^4;

% mean longitude of Jupiter (LJup)
arg{10} = 0.599547105074 + 52.993480508539*t + 0.000003897272*t.^2 + 0.000000000646*t.^3 -0.000000000091*t.^4;

% mean longitude of Saturn (LSat)
arg{11} = 0.874016284019 + 21.354296582042*t + 0.000009059625*t.^2 -0.000000000521*t.^3 -0.000000000170*t.^4;
