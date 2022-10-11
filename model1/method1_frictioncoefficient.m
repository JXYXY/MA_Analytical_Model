clc
clear
close all
% parameter
g = 9.81;
I_y = 0.008^3*0.015/12; % moment of area


rho = 2000; % line density of finger
Area = 0.008*0.015; % area of cross section 
L = 0.1; % length of finger
muMatrix = linspace(0,0.3,6);

H = 0.04; % X distance
E = exp(20*0.0235-0.6403)*10^6;
% E = 30*10^9;

for n = 1:1:6
    mu = muMatrix(n);
    % initialization
    F_N = 0; % grasping force
    delta_F_N = 0.01; % step size of grasping force
    T_d = 0; % grasping force
    delta_T_d = 0.001; % step size of grasping force
    e_H = 0.0001;
    X_p = 0;
    Z_p = L;

    flag = 0;

    while flag == 0

        A = -F_N/E/I_y;
        B = (mu*F_N + rho * Area * g * L)/E/I_y;
        C = - rho * Area * g /E/I_y;
        D = T_d/E/I_y;

        [x,z]= solFingerModel(A,B,C,D,L);
        X_p_new = x(end);
        Z_p = z(end);

        if X_p_new <= X_p && X_p_new > 0
            flag = 1;
        end
        X_p =X_p_new;

        if X_p >= H

            if abs(X_p - H) <= e_H
                T_d = T_d + delta_T_d;
            else
                F_N = F_N + delta_F_N;
            end
        else
            T_d = T_d + delta_T_d;

        end

    end

    F(n) = F_N;
    T(n)= T_d;
    
    plot(x,z)
    title("centerline position")
    xlabel X
    ylabel Y
    leg_str{n} = ['m=',num2str(mu)];
    hold on

end

legend(leg_str)
figure(2);
plot(muMatrix,F);
    title("gripping force")
    xlabel('friction coefficient')
    ylabel('gipping force')

    figure(3);
plot(muMatrix,T);
    title("moment")
    xlabel('friction coefficient')
    ylabel('moment')

