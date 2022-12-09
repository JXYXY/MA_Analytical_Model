function lambda = stretch(x,theta,dtheta,H,F_gr,mu,B,g_y,L)
syms l_iter
    eqn = l_iter - 1/l_iter^3 - 2*dtheta(10)^2*H^2/l_iter^9 == 0.5*F_gr(end)/mu/B/H*cos(theta(10))+(F_gr(end)+g_y*(L-x(10)))/mu/B/H*sin(theta(10));
    lambda_cal = vpasolve(eqn,l_iter,1e-4);
    lambda = double(lambda_cal(real(lambda_cal)>0&imag(lambda_cal)==0));
end
