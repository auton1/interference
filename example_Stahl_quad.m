function [nu_est, p_est, lk_max, covariance_matrix] = example_Stahl_quad(nu, p, cM_map_len, N_indv)
% Run a simulation study of the Stahl model (quad / phase unknown data).
%
% Usage: [nu_est, p_est, lk_max, covariance_matrix] = example_Stahl_quad(nu, p, cM_map_len, N_indv)
%
% Output:
% 	nu_est: estimate of interference parameter
%	p_est: estimate of escape parameter
%	lk_max: log likelihood of fitted model
%   covariance_matrix: covariance matrix estimated as inverse of fisher information matrix.
%
% Input:
%	nu: interference parameter to simulate
%   p: escape parameter to simulate
%	cM_map_len:	vector of chromosome map lengths in centiMorgans
%	N_indv: number of individuals to simulate
%

if (nargin < 2)
    error('Require a interference and escape parameter');
end

if (nargin < 3)
    % Simulate two chromosomes of lengths 200cM and 250cM
    cM_map_len = [200; 250];
end

if (nargin < 4)
    % Simulate 300 individuals
    N_indv = 300;
end

if ((nu < 0.1) | (nu > 50))
    warning('Function will not estimate nu outside range 0.1 to 50');
end

if ((p < 0) | (p > 0.5))
    error('p must be within range 0 to 0.5');
end

L = cM_map_len / 100;   % Convert to Morgans

opt = optimset ( 'Display', 'iter', 'TolX',1e-3);
events = simStahl_quad(N_indv, L, nu, p);

[res, lk_max] = fminsearchbnd(@(x)(-stahlLogLk_quad(events, L, x(1), x(2))), [1 eps], [0.1 eps], [50 0.5], opt);
lk_max = -lk_max;
nu_est = res(1);
p_est = res(2);

fisher_information = -hessian(@(x)(stahlLogLk_quad(events, L, x(1), x(2))), res);
covariance_matrix = inv(fisher_information); 


end
 
