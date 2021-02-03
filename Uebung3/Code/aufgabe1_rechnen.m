close all
clear all
clc
R = 6371000;
g = 2.72e-6;
gamma = 9.81;
h0 = 36.63;
tr0 = 357.9;
tr1 = 355.3;

T0 = R * g - 2*tr0 + 2 * gamma * h0;
dW = (g + T0/R)*R/2;

