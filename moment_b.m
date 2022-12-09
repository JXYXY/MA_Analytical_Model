function  [g_y,M_T] = moment_b(p,lambda,counter)

% geometrical parameters
% t1 = 0.005;
% a = 0.010; % h1/2
% b = 0.010; % breadth of air chamber 2b
% w1 = 0.020;
% d = 0.0025; % w2/2
% H = 0.010;
% N = 6;
% h_x = 0.0025;
% h_l = 0.020;
% L = 0.120;

t1 = 0.003;
a = 0.006; % h1/2
b = 0.006; % breadth of air chamber 2b
w1 = 0.015;
d = 0.0015; % w2/2
H = 0.006;
N = 7;
h_x = 0.0015;
h_l = 0.015;
L = 0.102;


% material properties
c1 = 0.0448e6;
c2 = 0.02154e6;
rho = 1000;
g = 9.81;

% line tensile in free inflation

theta0 = 1e-5;
fun_free_inflation = @(theta_mf)line_tensile(theta_mf,t1,c1,c2,p,a);
theta_mf = fsolve(fun_free_inflation,theta0);

% inclined contact
theta0 = [0.1 1];
fun_inclination = @(theta)inclination(theta,theta_mf,d,a,t1,c1,c2,p);
theta = fsolve(fun_inclination,theta0);
rho_star = a*tan(theta(1))/tan(theta_mf);
lambda_m = cos(theta(1)) + (d+2*a*sin(theta(1)))/2/(a+rho_star)*(theta(2)/sin(theta(2)/2)^2 - 2*cot(theta(2)/2));
R = (d+2*a*sin(theta(1)))/(2*sin(theta(2)/2)^2);
c_y = 2*a*lambda_m - R*(acos(1-d/R) + theta(2));
c_z = c_y*b*(a-rho_star)/a^2;
theta_s = theta(1);
% contact moment

% M_p = N * sqrt(2*h_x*M_m^2/(2*h_x+h_l));
% M_g = 2.3*rho*0.020^2*g*L^2/2;
V_chamber = (N*w1*2*a*2*b+(0.8*t1)^2*2*d*(N-1));
w0 = N*w1+(N-1)*2*d+2*N*t1;
d0 = 2*b+2*t1;
h0 = 2*a+t1;
V_gripper = w0*d0*(h0+2*t1)-V_chamber;


% M_g = (1+3*V_chamber/V_gripper)*rho*(0.030*0.035-0.020^2)*g*L^2/2;
% M_g = (1+0.656*2e-5*p)*rho*0.012^2*g*L^2/2;
    A_c = pi/4*c_y*c_z;
    F_p = p * A_c;
    ey = (H/2+a-rho_star)*cos(theta_s);
if c_y<0
    M_m = 0;
else

    M_m = ey * F_p;

end
M_p = N * M_m;    
if  c_y<0 && counter ~=0
V_chamber  = V_chamber *lambda;
else 
    V_chamber  = V_chamber *lambda_m;
end
g_y = (1+3*V_chamber/V_gripper)*g;
M_g = g_y*rho*2*a*2*b*L^2/2;
M_T =M_p + M_g;
end

