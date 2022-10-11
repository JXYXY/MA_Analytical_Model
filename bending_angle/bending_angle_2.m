clc
clear 
close all
%% unable to introduce external moment！！
% only with the wedged-shape chamber
% geometrical parameter 

W = 0.011;

b1 = 0.002;

l2 = 0.004;
L = 0.1;

H1 = 0.018;

N = 10;

F = 0.3;
M = F*L/2;
I_y = (W+2*b1)*(H1+2*b1)^3/12;

% material property ??

epsilon_x = @(sigma_x) 6.35e-7*sigma_x.^2 + 9.4e-3.*sigma_x - 0.017;


% loading without external moment
for n = 1:1:4
    thetaOfDegree = n*10;
    theta = thetaOfDegree*pi/180;

    H2 = H1 - l2/2/tan(theta);
    p = 1e4:1e4:5e4;

    sigma = p.*W*H1/(2*H1*b1+(W+2*b1)*b1) + M*H2/I_y;
    epsilon = epsilon_x(sigma);
    Dx = l2.*epsilon;
    beta = atan(Dx/2*H2);
    thetaOfChamber = beta/pi*180;
    thetaOfGripper = 2*N.*thetaOfChamber;
    R = L./beta;

    figure(1);
    plot(p/1000,thetaOfGripper);
    title("bending angle of gripper");
    xlabel("pressure [kPa]");
    ylabel("bending angle [degree]");
    leg_str_1{n} = ['theta=',num2str(thetaOfDegree)];
    hold on

    figure(2);
    plot(p/1000,R);
    title("bending radius of gripper");
    xlabel("pressure [kPa]");
    ylabel("bending radius [m]");
    leg_str_2{n} = ['theta=',num2str(thetaOfDegree)];
    hold on

    % figure(1);
    % plot(p/1000,thetaOfGripper)
    % title("bending angle of gripper");
    % xlabel("pressure [kPa]")
    % ylabel("bending angle [degree]")
end
legend(leg_str_1)
legend(leg_str_2)





