import numpy as np
import scipy.stats as stats
from functools import partial
from collections import defaultdict

# Simulate phase-known data under Stahl model.
# Following simStahl.c from xoi R package by Karl Broman.
#
# Usage: events = simStahl(n_sim, L, nu, p_escape)
#
# n_sim : number of simulations
# L :  np.array of map lengths in *Morgans*
# nu : interference parameter (float or np.array of length L)
# p_escape : probability that an event is not subject to 
#            interference (float or np.array of length L)
#
# events : n_sim by n_chr dictionary of event position lists, with one array
#          per simulation
#
def simStahl(n_sim, L, shape, p_escape=0.0):
	'''Simulate phase-known data under Stahl model.'''

	if type(L) is list: L = np.array(L)
	if type(L) is not np.ndarray: L = np.array([L]) 
	if type(shape) is not np.ndarray: shape = shape * np.ones(L.shape)
	if type(p_escape) is not np.ndarray: p_escape = p_escape * np.ones(L.shape)

	assert (n_sim > 0), 'n_sim must integer be greater than 0.'
	assert np.all(L > 0.0), 'L must be greater than 0.'
	assert np.all(shape > 0.0), 'shape must be greater than 0.'
	assert np.all(p_escape >= 0.0) & np.all(p_escape <= 1.0), 'p_escape must be between 0 and 1.'

	n_chr = len(L)
	if ((len(shape) == 1) & (len(L) > 1)):
		shape = shape*np.ones(np.shape(L))

	if ((len(p_escape) == 1) & (len(L) > 1)):
		p_escape = p_escape*np.ones(np.shape(L))

	scale = 1.0 / (2.0*shape*(1.0-p_escape))

	events = defaultdict(lambda : defaultdict(partial(np.ndarray, 0)))

	n_bins = 100000
	for c in range(n_chr):
		step = (1.0*L[c]) / n_bins
		bins = (np.array(range(n_bins))+0.5)*step
		startprob = np.cumsum(2.0 * (1.0 - p_escape[c]) * step * (1.0 - \
					stats.gamma.pdf(bins, shape[c], scale=scale[c])))

		for s in range(n_sim):
			events[s][c] = np.array([])
			u = np.random.rand()
			pos = 0.0
			if (u < startprob[-1]):
				# Choose a starting positon
				first_i = np.where(u < startprob)[0][0]
				for i in range(first_i, len(startprob)):
					if (u < startprob[i]):
						pos = (i-0.5)*step
						if (np.random.rand() < 0.5):
							events[s][c] = np.append(events[s][c], pos)
						break

				while (pos < L[c]):
					pos = pos + stats.gamma.rvs(shape[c], scale=scale[c])
					if ((np.random.rand() < 0.5) & (pos < L[c])):
						events[s][c] = np.append(events[s][c], pos)

			if (p_escape[c] > 0.0):
				# Add events that escape interference
				n_escape = stats.poisson.rvs(L[c] * p_escape[c])
				if (n_escape > 0):
					pos = L[c] * np.random.rand(n_escape)
					pos = np.concatenate((events[s][c], pos))
					pos = np.sort(pos)
					events[s][c] = pos
			events[s][c] = np.unique(events[s][c]) # Handle low prob of collision.
	return(events)

# Simulate phase unknown data from a quartet under the Stahl model. 
# Following simStahl.c from xoi R package by Karl Broman.
#
# Usage: events = simStahl_quad(n_sim, L, nu, p_escape);
#
# n_sim : number of simulations
# L : np.array of map lengths in *Morgans*
# nu : interference parameter (float or np.array of length L)
# p_escape : probability that an event is not subject to 
#			 interference (float or np.array of length L)
#
# events : n_sim by n_chr dictionary of event positions, with one array
#          per simulation
#
def simStahl_quad(n_sim, L, shape, p_escape=0.0):
	'''Simulate phase unknown data from a quartet under the Stahl model.'''
	if type(L) is not np.ndarray: L = np.array([L])
	if type(shape) is not np.ndarray: shape = shape * np.ones(L.shape)
	if type(p_escape) is not np.ndarray: p_escape = p_escape * np.ones(L.shape)

	assert (n_sim > 0), 'n_sim must integer be greater than 0.'
	assert np.all(L > 0), 'L must be greater than 0.'
	assert np.all(shape > 0), 'shape must be greater than 0.'
	assert np.all(p_escape >= 0) & np.all(p_escape <= 1), 'p_escape must be between 0 and 1.'

	events = simStahl(n_sim, L, shape, p_escape)
	events2 = simStahl(n_sim, L, shape, p_escape)

	n_chr = len(L)
	for c in range(n_chr):
		for s in range(n_sim):
			events[s][c] = np.sort(np.unique(np.concatenate((events[s][c], events2[s][c]))))
	return(events)
