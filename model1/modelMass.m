function [R] = modelMass(mMatrix,H,P,g,E,I_y,rho,Area,L)

T_d = 0.5*P*0.015*0.015^2;
for n = 1:numel(mMatrix)
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
        delta_theta = (theta - [0 theta(1:(end-1))]);
        X_p = x(end);


        Z_p = z(end);
        if X_p < max(x)
            Y= sprintf('pressure is too large, when m = %d p = %d H = %d',m,P,H);
            disp(Y)
            F(n) = 0;
            break;
        end
        if X_p < H
            X= sprintf('pressure is not enough, when m = %d p = %d H = %d',m,P,H);
            disp(X)
            F(n) = 0;
            break;

        end
        F_N = F_N - delta_F_N;
    end

    F(n) = F_N;
    pos = [x',z'];
    [K_L,R,K] = curvature(pos);

    subplot(1,4,1)
    plot(x*1000,z*1000)
    title("Gripper Position")
    xlabel('X [mm]')
    ylabel('Y [mm]')
    leg_str{n} = ['m=',num2str(m),'kg'];
    legend(leg_str)
    hold on

    subplot(1,4,2)
    numChamber = 1:numel(delta_theta);
    plot(numChamber,delta_theta)
    title("Bending angle of each chamber")
    xlabel('Chamber Number')
    ylabel('Bending Angle [deg]')
    leg_str{n} = ['m=',num2str(m)];
    txt = ['\leftarrow' num2str(theta(end)) 'deg in total'];
    text(numChamber(end-2),delta_theta(end-2),txt)
    legend(leg_str)
    hold on

    subplot(1,4,3)
    numChamber = 1:numel(R);
    plot(numChamber,R)
    title("Bending radius of each chamber")
    xlabel('Chamber Number')
    ylabel('Bending radius [deg]')
    leg_str{n} = ['m=',num2str(m)];
    legend(leg_str)
%     txt = ['\leftarrow' num2str(theta(end)) 'deg in total'];
%     text(numChamber(end-2),delta_theta(end-2),txt)
    hold on
end


subplot(1,4,4)

plot(mMatrix,-F)
title("Gripping Force")
xlabel('Mass [kg]')
ylabel('Gipping Force [N]')
hold on


end


