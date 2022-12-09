function  [x,theta,dtheta,F_gr,H_ext] = constitutive_equations_m(L,lambda,g,mu,B,H,C1)
x = linspace(0,L,50);
    solinit=bvpinit(x,[1e-5;1e-5;1e-5;1e-5]);
%     options = bvpset('RelTol',1e-6,'AbsTol',1e-10);
    sol = bvp5c(@(x,t)odefun(x,t,L,lambda,g,mu,B,H),@(xa,xb)bcfcn(xa,xb,C1),solinit);
    x = sol.x;
    theta = sol.y(1,:);
    dtheta = sol.y(2,:);
    F_gr = sol.y(end,:);
    H_ext = 0.5*F_gr;

end