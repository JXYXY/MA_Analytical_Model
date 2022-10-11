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

% muMatrix = 0.1:0.1:0.4;
mMatrix = [0.05 0.1 0.15];
% H = 0.03;
% P = 4e4;
HMatrix = [0.030 0.035 0.040];
PMatrix = [3.0e4 3.5e4 4.0e4];

E = exp(20*0.0235-0.6403)*10^6;

F = cell(numel(HMatrix)numel(PMatrix),numel(mMatrix));




for j = 1:numel(HMatrix)
    H = HMatrix(j);
    for k = 1:numel(PMatrix)
        P = PMatrix(k);
        T_d = 0.5*P*0.015*0.015^2;

        for n = 1:3
            m = mMatrix(n);
            F_g = m*g;
            %     mu = muMatrix(n);
            % initialization
            F_N = 0; % grasping force
            delta_F_N = 0.001; % step size of grasping force
            e_H = 0.0001;
            X_p = 0;
            Z_p = L;
            while abs(X_p - H) > e_H

                A = -F_N/E/I_y;
                B = (F_g/2 + rho * Area * g * L)/E/I_y;
                C = - rho * Area * g /E/I_y;
                D = T_d/E/I_y;

                [x,z,theta]= solFingerModel(A,B,C,D,L);
                delta_theta = (theta - [0 theta(1:(end-1))])*180/pi;
                X_p = x(end);


                Z_p = z(end);
                if X_p < max(x)
                    Y= sprintf('pressure is too large, when m = %d p = %d H = %d',m,P,H);
                    disp(Y)
                    F{} = 0;
                    break;
                end
                if X_p < H
                    X= sprintf('pressure is not enough, when m = %d p = %d H = %d',m,P,H);
                    disp(X)
                    F{j,k,n} = 0;
                    break;

                end
                F_N = F_N - delta_F_N;
            end

            F{j,k,n} = F_N;
            pos = [x',z'];
            [K_L,R,K] = curvature(pos);


            figure(j)
            subplot(1,4,1)
            plot(x,z)
            title("centerline position")
            xlabel X
            ylabel Y
            leg_str{n} = ['m=',num2str(m)];
            hold on

            subplot(1,4,2)

            plot(x,theta*180/pi,x,delta_theta)
            title("bending angle")
            xlabel X
            ylabel Y
            leg_str{n} = ['m=',num2str(m)];
            hold on


        end
    end
end

legend(leg_str)
subplot(1,4,3)

plot(mMatrix,-F)
title("gripping force")
xlabel('mass')
ylabel('gipping force')
hold on

subplot(1,4,4)
plot(mMatrix,T);
title("moment")
xlabel('mass')
ylabel('moment')
hold on




