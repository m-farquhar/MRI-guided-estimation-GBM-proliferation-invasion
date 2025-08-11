function C = analytic_sol_gauss_C_nondim(R, tau, alpha)
% Computes the nondimensional concentration C at nondimensional radius R,
% nondimensional time tau, with nondimensional parameter alpha.

C = exp(tau - R.^2./(4*alpha*tau+1))./(4*alpha*tau+1).^(3/2);