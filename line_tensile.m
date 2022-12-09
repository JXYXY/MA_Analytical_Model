function T = line_tensile(theta_mf,t1,c1,c2,p,a)
lambda = theta_mf/sin(theta_mf);
T = 2*t1*(c1+c2)*(lambda - 1/lambda^3)-p*a/sin(theta_mf);
% T = 2*t1*(c1+c2)*(1 - 1/lambda^6)-p*a/sin(theta_mf);
end