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
lambda_z = (9 + 11/60) / 180 * pi;
phi_z = (48 + 46/60) /180 * pi;

% P(6,4)
theta = linspace(0,pi,181);
P_all = Normalized_Lengendre(6,theta);
P_64 = P_all{7,5};

% Y
lambda = linspace(-pi,pi,361);
% to be continued


