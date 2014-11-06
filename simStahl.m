function events = simStahl(n_sim, L, shape, p_escape)
% Simulate phase-known data under Stahl model. Following simStahl.c from 
% xoi R package by Karl Broman.
%
% Usage: events = simStahl(n_sim, L, nu, p_escape);
%
% n_sim : number of simulations
% L : vector of map lengths in *Morgans*
% nu : interference parameter
% p_escape : probability that an event is not subject to interference
%
% events : n_sim by n_chr cell array of event positions, with one cell per 
%            simulation
%

if (nargin < 4)
    p_escape = 0;   % Simulate under Gamma Model without escape.
end

assert((n_sim > 0), 'n_sim must integer be greater than 0.');
assert(all(L > 0), 'L must be greater than 0.');
assert(all(shape > 0), 'shape must be greater than 0.');
assert(all(p_escape >= 0) && all(p_escape <= 1), 'p_escape must be between 0 and 1.');

N_chr = length(L);

if ((length(shape) == 1) & (length(L) > 1))
    shape = shape*ones(size(L));
end

if ((length(p_escape) == 1) & (length(L) > 1))
    p_escape = p_escape*ones(size(L));
end

scale = 1./(2.*shape.*(1-p_escape));

events = cell(n_sim, N_chr);

N_bins = 10000;
for c=1:N_chr
    step = L(c) / N_bins;
    bins = ((0:(N_bins-1))+0.5)*step;

    startprob = cumsum(2.0 * (1 - p_escape(c)) * step * (1-gamcdf(bins, shape(c), scale(c))));
    
    chr_events = events(:,c);

    parfor s=1:n_sim
        u = rand;
        pos = 0;
        if (u < startprob(end))
            %chr_events{s} = [];
            % Choose a starting position
            first_i = find(u<startprob, 1, 'first');
            for i=first_i:length(startprob)
                if (u < startprob(i))
                    pos = (i-0.5)*step;
                    if (rand < 0.5)
                        chr_events{s}(end+1) = pos;
                    end
                    break;
                end
            end

            while(pos < L(c))
                pos = pos + gamrnd(shape(c), scale(c));
                if ((rand < 0.5) && (pos < L(c)))
                    chr_events{s}(end+1) = pos;
                end
            end
        end

        if (p_escape > 0)
            % Add events that escape interference
            n_escape = poissrnd(L(c) * p_escape(c));
            if (n_escape > 0)
                pos = L(c)*rand(1, n_escape);
                chr_events{s} = sort([chr_events{s} pos]);
            end
        end
    end
    
    events(:,c) = chr_events;
end

