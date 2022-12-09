function dtdx = odefun(x,t,L,lambda,g,mu,B,H)
% (33)(43)(47) ODE fUNCTION
rho = 1000;
g_y =g*rho*B*H;
c = (mu*B*H^3)/3/lambda^8;

dtdx = zeros(4,1);
dtdx = [
    t(2);
    -t(4)*0.5/c * sin(t(1)) +  t(4)/c*cos(t(1)) - g_y/c*(L-x)*cos(t(1));
    lambda*sin(t(1));
    0
    ];
end