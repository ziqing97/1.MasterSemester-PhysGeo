%%
clc
% Nicholas Schneider & Ziqing Yu
% 25/11/2020

% Aufgabe 2b test
% f2
[lambda,theta] = meshgrid(linspace(-pi,pi,200), linspace(0,pi,100));

t21 = sin(theta).^3 .* sin(lambda);
t22 = 0.75 * sin(theta) .* sin(lambda) - 0.25 * sin(3*theta) .* sin(lambda);
d2 = t21 - t22;

% f3
t31 = cosh(cos(theta));
t32 = sinh(1) + (20 * sinh(1) - 15 * cosh(1)) * 0.25 * (1+3*cos(2 * theta)) + (453 * sinh(1) - 345 * cosh(1)) * 3 ...
    *(35/8*cos(theta).^4 - 15/4*cos(theta).^2 + 3/8);
d3 = t31 - t32;

