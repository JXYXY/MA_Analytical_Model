clc
clear 
close all

% geometrical parameter !! may lead to failure 
% V and I_y differentassumption
H = 0.018;
H1 = 0.011;
H2 = 0.002;
H3 = 0.005;
Hc = H1*0.5;  % change Hc for SPNA and FPNA

W = 0.015;
W1 = 0.002;
Wc = W - 2*W1;
L20 = 0.003;

N = 10;

% H = 0.020;
% H1 = 0.010;
% H2 = 0.002;
% H3 = 0.008;
% Hc = H1/2;
% W = 0.014;
% W1 = 0.002;
% Wc = W - 2*W1;
% L20 = 0.002;
% N = 14;


I = W*H3^3/12 - Wc*H1^3/12;


% material property
% E = exp(20*0.0235-0.6403)*10^6;
E = exp(20*0.0235-0.6403)*10^6;

% loading without external moment
p = 0:1e4:3e4;

k = ( E*I + sqrt((E*I)^2 - 1/2*pi*E*I*Wc*H1*Hc^2.*p)) / (E*I*Hc); % sign? 
R = 1./k;
thetaOfChamber =  L20./(R - Hc)/pi*180;
thetaOfGripper = N.*thetaOfChamber;
plot(p,thetaOfGripper)


