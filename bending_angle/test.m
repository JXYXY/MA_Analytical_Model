clc
clear 
close all

% geometrical parameter !! may lead to failure 
% V and I_y different assumption
H = 0.018;
H1 = 0.011;
H2 = 0.002;
H3 = H - H1 - H2;
Hc = H1*0.8;  

W = 0.015;
W1 = 0.002;
Wc = W - 2*W1;
L20 = 0.004;

N = 10;

I = W*H^3/12 - Wc*H1^3/12;

% material property
E = exp(20*0.0235-0.6403)*10^6;

% loading with external moment
p = 4e4;
F = 0.55;

syms k
L2(k) = L20/(1-Hc*k);
L = 0:0.008:10*0.008;

RSol = zeros(11,2);
ThetaSol = zeros(11,2);

for n = 2:1:11
l1 = L(n-1);
l2 = L(n);
sigma(k) = 1/2*E*I*k^2*L2 - p*pi/4*H1*L2*Wc- 1/2*(l2^2-l1^2)*F*k;
eqn = diff(sigma(k)) ==0;
kSol = vpasolve(eqn);
R = 1./kSol;
RSol(n,:) = R;
theta = round(L2(kSol)./R*180/pi,3);
ThetaSol(n,:) = theta;
end

sol = sum(ThetaSol);
num = 1:11

bendingRadiusChamber=RSol(:,1)
figure(1);
plot(num',bendingRadiusChamber)


bendingAngleChamber=ThetaSol(:,1)
figure(2);
plot(num',bendingAngleChamber)

index=find(sol>0);
bendingAngleGripper=sol(:,1)



