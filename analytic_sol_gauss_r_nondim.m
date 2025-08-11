function R = analytic_sol_gauss_r_nondim(alpha, tau, beta)
% Computes the nondimensional detectable radius R for nondimensional parameter 
% alpha, at nondimensional time tau, with nondimensional detection threshold beta. 
% Anything equal to 0 (i.e. invisible), is set to a value of -1 so that it does 
% not show on the plots.

a = (tau - log((4*alpha*tau+1).^(3/2)*beta)).*(4*alpha*tau+1);

a(a<0) = 0;

R = sqrt(a);

R(R==0) = -1;
