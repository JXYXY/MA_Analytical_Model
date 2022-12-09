function eqns = inclination(theta,theta_mf,d,a,t1,c1,c2,p)

% theta_s = theta(1)
% theta_m2 = theta(2)

rho_star = a*tan(theta(1))/tan(theta_mf);
lambda = cos(theta(1)) + (d+2*a*sin(theta(1)))/2/(a+rho_star)*(theta(2)/sin(theta(2)/2)^2 - 2*cot(theta(2)/2));
R = (d+2*a*sin(theta(1)))/(2*sin(theta(2)/2)^2);
T = 2*t1*(c1+c2)*(lambda-1/lambda^3);
% T = 2*t1*(c1+c2)*(1-1/lambda^6);
eqns = [
    csc(theta_mf) - cos(theta(1))*cot(theta_mf) - sin(theta(1)) - d/a
    T/R - p];

end