clc
clear
close all
% parameter
g = 9.81;
HChamber = 0.015;
WChamber = 0.015;
LChamber = 0.008;
I_y = WChamber*HChamber^3/12; % moment of area
rho = 2000;% density of finger
Area = WChamber*HChamber;% area of cross section
L = 0.1; % length of finger

E = exp(20*0.0235-0.6403)*10^6;

model = 3;

switch model
    case 1
        mMatrix = [0.05 0.1 0.15];
        H = 0.03;
        P = 4e4;
        modelMass(mMatrix,H,P,g,E,I_y,rho,Area,L);
    case 2
        m = 0.1;
        HMatrix = [0.030 0.035 0.040];
        P = 4e4;
        modelSize(HMatrix,m,P,g,E,I_y,rho,Area,L);
    case 3
        m = 0.1;
        H = 0.03;
        PMatrix = [3.0e4 3.5e4 4.0e4];
        modelPressure(PMatrix,H,m,g,E,I_y,rho,Area,L);
end



