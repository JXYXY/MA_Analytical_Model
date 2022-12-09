clc
clear
close all

% L = 0.145; % length of actuator
% B = 0.010;
% H = 0.020;
L = 0.102;
B = 0.012;
H = 0.012;

g = 9.81;

rho = 1000;
c1 = 0.0448e6;
c2 = 0.02154e6;
mu = 2*(c1+c2); % shear modulus
g_y = g*rho*B*H;
% p_matrix = 5.6e4:2e3:8e4;
p_matrix = 5.2e4:2e3:8e4;
F_gr_matrix = zeros(1,length(p_matrix));
lambda_matrix = zeros(1,length(p_matrix));

for i = 1:length(p_matrix)
    counter = 0;
    lambda = 1;

    p = p_matrix(i);
    [g_y,M_T] = moment_b(p,lambda,counter);
    C1 = 3*M_T*lambda^7/(mu*B*H^3);

    % initial value

    [x,theta,dtheta,F_gr,H_ext] = constitutive_equations_m(L,lambda,g_y,mu,B,H,C1);
    theta_calo = theta;
    [g_y,M_T] = moment_b(p,lambda,counter);
    M_T= M_T-F_gr(end)*L;


    %% eqn 33
    lambda = stretch(x,theta,dtheta,H,F_gr,mu,B,g_y,L);
    theta_caln = zeros(size(theta_calo));

    
    while abs(theta_caln(end)-theta_calo(end)) > 0.0001
        theta_calo = theta_caln;
        [x,theta,dtheta,F_gr,H_ext] = constitutive_equations_m(L,lambda,g,mu,B,H,C1);
        theta_caln = theta;

        [g_y,M_T] = moment_b(p,lambda,counter) ;
        M_T= M_T-F_gr(end)*L;
        lambda = stretch(x,theta,dtheta,H,F_gr,mu,B,g_y,L);
        counter =counter+1;
    end
    F_gr_matrix(i) = F_gr(end);
    lambda_matrix(i) = lambda;

end
% data = xlsread( 'simulation_6chamber_10_20_20.xlsx' ) ;
data = xlsread( 'simulation_7chamber_6_12_12.xlsx' ) ;
plot(data(:,1),data(:,2))
hold on
plot(p_matrix,F_gr_matrix,'o')
xlabel('Pressure[Pa]')
ylabel('Gripping Force[N]')
legend('simulation','analytical model')
