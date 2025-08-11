function r = analytic_sol_gauss_r_dim(t, D, rho,b,cstar,c0)
% Computes the dimensional detectable radius r for dimensional parameters 
% D, rho, b, at dimensional time t, with detection threshold cstar. 
% Anything equal to 0 (i.e. invisible), is set to a value of -0.001 so that it does 
% not show on the plots.

a = (rho.*t - log((cstar*(b.^2 + 4*D.*t).^(3/2))./(c0*b.^3)));

r = sqrt(a.*(b.^2 + 4*D.*t));

r(a<1e-4) = -0.001;

