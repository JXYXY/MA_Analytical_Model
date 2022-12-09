clc
clear 
close all
% INPUT: ME,DE,AE,BOUNDARY CONDITIONS
% OUTPUT: THETA,FORCE

% Parameters
h_x = 2.5e-3;
h_L = 20e-3;
t_m = 5e-3;
a = 8e-3;
b = 20e-3; %??
H = 9e-3;
p = 7e4;
mu = 1e5;
L = 145e-3;
rho = 1000;
A = h_L*b;
g = 9.81;

F = 0;

% inclination fixed
theta_s = 0;
% M_T
d = h_x;
theta_c = slope(d,a,theta_s);
rho_star = a*theta_s/theta_c;
[lambda_m, theta_m2] = inflation(p,d,mu,t_m,a,theta_s,rho_star);
R = (d+2*a*sin(theta_s))/2/sin(theta_m2/2)^2;
c2 = (a + rho_star)*cos(theta_s) - (d + 2*a*sin(theta_s))*cot(theta_m2/2);
c1 = 2*a*lambda_m - R*(acos(1-d/R) + theta_m2)-c2;
c_y = c1 + c2;
c_z = c_y * b * (a-rho_star)/a^2;

A_c = pi/4 * c_y * c_z;
e_y = (H/2 + a - rho_star)*cos(theta_s);
F_p = p * A_c;

M_m = F_p * e_y;
M_P = (2*h_x*M_m^2/(2*h_x+h_L))^0.5;
M_g = rho*A*g*L^2/2;
M_F = F*L;
M_T = 10*M_P + M_g - M_F;


lambda = 1;
H_ext = 0;
F_gr = 0;
[theta,d_theta] = diffeqn(L, lambda,H_ext,F_gr,g,mu,b,H)

