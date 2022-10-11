function [x,z,theta]= solFingerModel(A,B,C,D,L)

% % numerically solving finger model
% 
% smesh = linspace(0,L,15);
% solinit = bvpinit(smesh, @guess);
% sol = bvp4c(@bvpfcn, @bcfcn, solinit);
% % x = sol.x;
% % z = sol.y(1,:);
% % delta_s = [sol.x(1) diff(sol.x)];
% % x = cumsum(delta_s.* sin(sol.y(1,:)));
% % z = cumsum(delta_s.* cos(sol.y(1,:)));

% numerically solving finger model
opts = bvpset('FJacobian',@jac,'RelTol',0.1,'AbsTol',0.1);
smesh = linspace(0,L,10);
solinit = bvpinit(smesh, [2*L;2*L]);
sol = bvp4c(@bvpfcn, @bcfcn, solinit,opts);
% x = sol.x;
% z = sol.y(1,:);
delta_s = [sol.x(1) diff(sol.x)];
x = cumsum(delta_s.* sin(sol.y(1,:)));
z = cumsum(delta_s.* cos(sol.y(1,:)));
theta = sol.y(1,:)*180/pi;

function dthetads = bvpfcn(s,theta)
dthetads = [ theta(2);...
    A * cos(theta(1)) + ( B + C * s) * sin(theta(1))];
end
%--------------------------------
function res = bcfcn(thetaa,thetab) % boundary conditions
res = [
       thetaa(1)
       thetab(2)-D];
end
%--------------------------------
% function g = guess(s) % initial guess for y and y'
% g = [
%     sin(s)
%     cos(s)
%     ];
% %     -L/pi*sin(pi/L*s)
% %     -cos(pi/L*s)];
% end
%--------------------------------
function dfdtheta = jac(s,thetha) % initial guess for y and y'
dfdtheta = [0 1
    -A*sin(thetha(1))+(B+C*s)*cos(thetha(1)) 0
    ];

end


end
