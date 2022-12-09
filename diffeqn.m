function [theta,d_theta,y] = diffeqn(L,lambda,H_ext,F_gr,g,mu,B,H)
    x = linspace(0,L,5);
    solinit = bvpinit(x, [0.001;0.001;0.001]);
    M_T = 2.0431E4;
    sol = bvp4c(@odefun, @bcfcn,solinit);
    theta = sol.y(1,:);
    d_theta = sol.y(2,:);
    y = sol.y(3,:);





end