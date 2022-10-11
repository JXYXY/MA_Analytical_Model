% clc
% clear
% close all

syms e11 e22 e33 k1 k2 q x3

% geometrical parameter 
t1 = 1e-3;
t2 = 1e-3; 
t3 = 3e-3; 

h1 = 3.5e-3;
h2 = 2e-3;
h3 = 13e-3;
h4 = 2e-3;
h5 = 2e-3;

w = 15e-3; % chamber width 

theta = pi/2; % chamber angle
phi = 0; % twisting angle
L = 0.1; % gripper length
n = 10; % number of chambers
I_Y = 0.03^3*0.05/12; %moment of area

% material properties
E = exp(20*0.0235-0.6403)*10^6; % Young's modulus
nu = 0.5; % Poisson ratios
c = 1;  % effective parameter

% loading
P = 200000; % pressure [pa]

% strain energy density function ws

sigma = diag([e11 + x3*k1 , e22 + x3*k2 , e33 + x3*q]); % strain
J1 = trace(sigma); % first strain invariant
J2 = trace(sigma^2); % second strain invariant
ws = 1/2*E/(1+nu)*J2 + E*nu*(1-2*nu)/(1+nu)*J1^2; 

% elastic energy W_m

wm1 = 2*w*(t1+t2+t3)*int(ws,x3,0,h1);
wm2 = 2*w*(t1+t2)*int(ws,x3,h1+h2,h1+h2+h5);
wm3 = 2*w*(t2)*int(ws,x3,h1+h2+h5,h1+h2+h3);
wm4 = 2*w*(t2+t3)*int(ws,x3,h1+h2+h3,h1+h2+h3+h4);
W_m = wm1+wm2+wm3+wm4;

% potential energy of force W_f
delta_W_f = c * 2* P *  w * t3/sin(theta)*((sigma(1,1)*sin(theta-phi)^2+sigma(2,2)*cos(theta-phi)^2));
W_f = int(delta_W_f,x3,h1+h2,h1+h2+h3);

% potential energy
potentialEng = n*(W_m - W_f);

% minimum potential energy

f = [diff(potentialEng, e11)
    diff(potentialEng, e22)
    diff(potentialEng, e33)
    diff(potentialEng, k1)
    diff(potentialEng, k2)
    diff(potentialEng, q)]==0;

x=solve(f,[e11,e22,e33,k1,k2,q]);

% centerline 
alpha = sqrt(x.k1^2*cos(phi)^2+x.k2^2+sin(phi)^2);
beta = x.k1*cos(phi)^2+x.k2^2+sin(phi)^2;
tao = (x.k2-x.k1)*sin(phi)*cos(phi);
m= length(0:0.01:L);
X1 = zeros(1,m);
X2 = zeros(1,m);
X3 = zeros(1,m);
for s = 0:0.01:L
    X1(s/0.01+1) = s * cos(phi)^2 - x.k1*s*beta*cos(phi)^2/alpha^ 2 +...
        x.k1*beta*cos(phi)^2*sin(s*alpha)/alpha^3 + s*sin(phi)^2 + ...
        - x.k2*s*beta*sin(phi)^2/alpha^ 2 + x.k2*beta*sin(phi)^2*sin(s*alpha)/alpha^3;
    X2(s/0.01+1) = tao *(x.k1 + x.k2 + (x.k1-x.k2)*cos(2*phi))*(s*alpha - sin(s*alpha)) /2/ alpha^3;
    X3(s/0.01+1) = (cos(s*alpha) -1)*(x.k1+x.k2+(x.k1-x.k2)*cos(2*phi))/2/ alpha^2;
end
figure;
plot3(X1,X2,X3)
title('Centerline Position')
xlabel('X')
ylabel('Y')
zlabel('Z')



% euler-bernoulli beam equation to calculate gripping force
X = [X1',X2',X3'];
[K_L,R,K] = curvature(X);
figure;
plot(K_L,1./R)
xlabel K_L
ylabel R
title('Curvature vs. cumulative curve length')

M = E*I_Y./R;
F_gripping= M/L;
ind = find(~isnan(F_gripping));
F_gripping = F_gripping(ind);

if sum(abs(diff(F_gripping))) < 0.001
    sprintf("gripping force is %.3f N",F_gripping(end))
else
    disp("Error in calculating gripping force!");
end


% end effector displacement
delta_d = norm([X1(end) X2(end) X3(end)] - [L 0 0])';

sprintf("end effector displacement is %4f m",delta_d)

    
