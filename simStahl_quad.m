function events = simStahl_quad(n_sim, L, shape, p_escape)
% Simulate quartet data under Stahl model. Following simStahl.c from xoi 
% package by Karl Broman.
%
% Usage: events = simStahl_quad(n_sim, L, nu, p_escape);
%
% n_sim : number of simulations
% L : vector of map lengths in *Morgans*
% nu : interference parameter
% p_escape : probability that an event is not subject to interference
%
% events : n_sim by n_chr cell array of event positions, with one cell 
%          per simulation
%

if (nargin < 4)
    p_escape = 0;   % Simulate under Gamma Model without escape.
end

assert((n_sim > 0), 'n_sim must integer be greater than 0.');
assert(all(L > 0), 'L must be greater than 0.');
assert(all(shape > 0), 'shape must be greater than 0.');
assert(all(p_escape >= 0) && all(p_escape <= 1), 'p_escape must be between 0 and 1.');

events = simStahl(n_sim, L, shape, p_escape);
events2 = simStahl(n_sim, L, shape, p_escape);

N_chr = length(L);
for c=1:N_chr
    for s=1:n_sim
        events{s,c} = sort([events{s,c} events2{s,c}]);
    end
end