function c = analytic_sol_gauss_C_dim(r, t, D, rho,b,c0)
% Computes the dimensional concentration c at radius r, time t, with parameters
% diffusion D, proliferation rho, initial spread b and initial
% concentration c0.

c=c0*exp(rho*t-r.^2./(b^2 + 4*D*t))./((b^2 + 4*D*t).^(3/2)) * b^3;


