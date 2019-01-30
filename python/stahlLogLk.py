import numpy as np
import sys
import scipy.stats as stats
import scipy.integrate as integrate
import pdb

def fstar1(x, p):
	if type(x) is not np.ndarray: x = np.array([[x]])
	if x.ndim == 1:
		x = np.array([x])
	k_max = 25
	k_div = np.power(0.5, range(1,k_max+1)).T
	a = np.dstack((x,) * k_max)
	b = np.tile(np.reshape(np.array(range(1, k_max+1)), (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	c = np.tile(np.reshape(k_div, (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	tmp = stats.gamma.pdf(a, b, scale=0.5/p) * c
	return np.sum(tmp, axis=2)

def Fstar1(x,p):
	if type(x) is not np.ndarray: x = np.array([[x]])
	if x.ndim == 1:
		x = np.array([x])
	k_max = 25
	k_div = np.power(0.5, range(1,k_max+1)).T
	a = np.dstack((x,) * k_max)
	b = np.tile(np.reshape(np.array(range(1, k_max+1)), (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	c = np.tile(np.reshape(k_div, (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	tmp = stats.gamma.cdf(a, b, scale=0.5/p) * c
	return np.sum(tmp, axis=2)	

def gstar1(x,p):
	return p * (1.0 - Fstar1(x,p))

def Gstar1(x,p):
	out = integrate.quad(lambda X:gstar1(X,p), 0, x)
	return out[0]

def fstar2(x,nu,p):
	if type(x) is not np.ndarray: x = np.array([[x]])
	if x.ndim == 1:
		x = np.array([x])
	k_max = 25
	k_div = np.power(0.5, range(1,k_max+1)).T
	a = np.dstack((x,) * k_max)
	b = np.tile(np.reshape(np.array(range(1, k_max+1))*nu, (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	c = np.tile(np.reshape(k_div, (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	tmp = stats.gamma.pdf(a, b, scale=0.5/((1-p)*nu)) * c
	return np.sum(tmp, axis=2)

def Fstar2(x,nu,p):
	if type(x) is not np.ndarray: x = np.array([[x]])
	if x.ndim == 1:
		x = np.array([x])
	k_max = 25
	k_div = np.power(0.5, range(1,k_max+1)).T
	a = np.dstack((x,) * k_max)
	b = np.tile(np.reshape(np.array(range(1, k_max+1))*nu, (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	c = np.tile(np.reshape(k_div, (1,1,k_max)).T, (x.shape[1], x.shape[0])).T
	tmp = stats.gamma.cdf(a, b, scale=0.5/((1-p)*nu)) * c
	return np.sum(tmp, axis=2)

def gstar2(x,nu,p):
	return (1.0-p)*(1.0 - Fstar2(x,nu,p))

def Gstar2(x,nu,p):
	out = integrate.quad(lambda X:gstar2(X,nu,p), 0, x)
	return out[0]

def dec2base3(d, width=0):
	if type(d) is not np.ndarray: d = np.array([d])
	n = int(np.max((1, np.rint(np.log2(np.max(d)+1.0)/np.log2(3.0)))))
	while np.any(np.power(3.0, n) <= d):
		n = n + 1
	if (width > 0):
		assert n <= width, 'width too short for required number'
	s = np.zeros(n)
	s[n-1] = np.fmod(d,3)
	while (np.any(d) & (n > 1)):
		n -= 1
		d = np.floor(d/3.0)
		s[n-1] = np.fmod(d,3.0)
	if (width > len(s)):
		s = np.pad(s, (width-len(s), 0), 'constant')
	return s

# Calculate log likelihoods according to the Stahl model.
#
# Events is a 2d dictionary of np.arrays.
# Each entry is a vector of crossover positions (in Morgans) from a given 
# individual on a given chromosome.
#
# L is the chromosome length (in Morgans).
# nu is the interference parameter.
# p is the escape probability.
#
# Note this function is likely to be very difficult to make sense of 
# as it has been highly optimized.
#
def stahlLogLk(events, L, nu, p):
	'''Calculate the log likelihoods according to the Stahl model.'''
	if type(L) is not np.ndarray: L = np.array([L])

	if (nu <= 0):
		return -99999999.9+nu
	if (nu > 100):
		return -99999999.9-nu
	if ((p < 0) | (p > 1)):
		return -99999999.9+p

	p = max(p, sys.float_info.epsilon)

	n_indv = len(events.keys())
	n_chr = len(events[0])
	assert n_chr == len(L), 'events does not match L in size'
	n_tot = n_indv * n_chr
	lk = 0.0

	log_1_minus_Gstar = np.zeros(n_tot)
	log_1_minus_Gstar_escape = np.zeros(n_tot)
	max_map = np.zeros(n_tot)
	all_events = [events[x][y] for x in events.keys() for y in events[x].keys()]

	tmp1 = np.zeros(n_chr)
	tmp2 = np.zeros(n_chr)
	for c in range(n_chr):
		tmp1[c] = np.log(1.0 - Gstar2(L[c], nu, p))
		tmp2[c] = np.log(1.0 - Gstar1(L[c], p))

	for c in range(n_chr):
		idx = np.array(range((c*n_indv),((c*n_indv)+n_indv)))
		log_1_minus_Gstar[idx] = tmp1[c]
		log_1_minus_Gstar_escape[idx] = tmp2[c]

	n_partitions = np.array([np.power(2, len(x)) for x in all_events])
	uniq_partitions = np.unique(n_partitions)
	n_partitions_reshape = np.reshape(n_partitions, (n_indv, n_chr))

	idx = n_partitions == 1
	lk = lk + np.sum(log_1_minus_Gstar[idx] + log_1_minus_Gstar_escape[idx])

	for u in range(1, len(uniq_partitions)): # Could parallelize this loop
		n_part = uniq_partitions[u]	

		width = len(np.binary_repr(n_part-1))
		partitions = [np.binary_repr(x, width=width) for x in range(n_part)]

		idx_shape = np.where(n_partitions_reshape == n_part)
		idx = np.where(n_partitions == n_part)
		all_pos = np.array([events[idx_shape[0][x]][idx_shape[1][x]] for x in range(len(idx_shape[0]))])
		lka = np.tile(log_1_minus_Gstar[idx], (n_part, 1))
		lkc = np.tile(log_1_minus_Gstar_escape[idx], (n_part, 1))
		max_map_idx = L[idx_shape[1]]
		
		I0 = np.array([list(x) for x in partitions]) == '0'
		I1 = I0 == False
		# Only consider unique partitions
		uniq_I0 = I0
		sum_sites0 = np.sum(uniq_I0,1)
		n_sites0 = uniq_I0.shape[1]
		uniq_I1 = I1
		sum_sites2 = np.sum(uniq_I1,1)

		#Loop over partitions with the same number of sites
		for n in range(1,n_sites0+1):
			J0 = sum_sites0 == n
			uniq_I00 = uniq_I0[J0,:]
			# Append boundaries
			uniq_I000 = np.concatenate((np.ones(uniq_I00.shape[0], dtype=bool)[:,np.newaxis], uniq_I00, np.ones(uniq_I00.shape[0], dtype=bool)[:,np.newaxis]), axis=1)
			all_pos0 = np.tile(np.concatenate((np.zeros(all_pos.shape[0])[:,np.newaxis], all_pos, max_map_idx[:,np.newaxis]), axis=1), (uniq_I000.shape[0], 1))

			uniq_I0000 = np.reshape(np.tile(uniq_I000, (1, all_pos.shape[0])), all_pos0.shape)
			all_pos00 = np.ma.MaskedArray(all_pos0, ~uniq_I0000)
			x0 = np.diff([x.compressed() for x in all_pos00], axis=1)
			
			y0 = np.log(gstar2(x0[:,0], nu, p) * (1.0 - Fstar2(x0[:,-1], nu, p)))
			lka[J0,:] = np.reshape(y0, lka[J0,:].shape)
			
			if (n > 1): # Check
				y1 = np.log(np.prod(fstar2(x0[:,1:-1], nu, p), axis=1))[:, np.newaxis]
				lka[J0,:] = lka[J0,:] + np.reshape(y1, lka[J0,:].shape)

			J1 = sum_sites2 == n
			uniq_I11 = uniq_I1[J1,:]
			# Append boundaries
			uniq_I111 = np.concatenate((np.ones(uniq_I11.shape[0], dtype=bool)[:,np.newaxis], uniq_I11, np.ones(uniq_I11.shape[0], dtype=bool)[:,np.newaxis]), axis=1)
			all_pos1 = np.tile(np.concatenate((np.zeros(all_pos.shape[0])[:,np.newaxis], all_pos, max_map_idx[:,np.newaxis]), axis=1), (uniq_I111.shape[0], 1))

			uniq_I1111 = np.reshape(np.tile(uniq_I111, (1, all_pos.shape[0])), all_pos1.shape)
			all_pos11 = np.ma.MaskedArray(all_pos1, ~uniq_I1111)
			x2 = np.diff([x.compressed() for x in all_pos11], axis=1)
			
			y2 = np.log(gstar1(x2[:,0], p) * (1.0 - Fstar1(x2[:,-1], p)))
			lkc[J1,:] = np.reshape(y2, lkc[J1,:].shape)
			
			if (n > 1): # Check
				y3 = np.log(np.prod(fstar1(x0[:,1:-1], p), axis=1))[:, np.newaxis]
				lkc[J1,:] = lkc[J1,:] + np.reshape(y3, lkc[J1,:].shape)
		
		tmp = lka + lkc
		ma = np.max(tmp, axis=0)
		tmp1 = ma + np.log(np.sum(np.exp(tmp - ma), axis=0))
		lk = lk + np.sum(tmp1)
	lk = lk + np.log(1.0 - p)
	print ("nu = %s p = %s lk = %s" % (nu, p, lk))
	return lk

# Calculate log likelihoods according to the Stahl model, but
# assuming that the events are observed in a quartet (i.e. are phase unknown).
#
# Events is a 2d dictionary of np.arrays.
# Each entry is a vector of crossover positions (in Morgans) from a given 
# individual on a given chromosome.
#
# L is the chromosome length (in Morgans).
# nu is the interference parameter.
# p is the escape probability.
#
# Note this function is likely to be very difficult to make sense of 
# as it has been highly optimized.
#
def stahlLogLk_quad(events, L, nu, p):
	'''Calculate the log likelihoods according to the Stahl model for quartet events.'''
	if type(L) is not np.ndarray: L = np.array([L])

	if (nu <= 0):
		return -99999999.9+nu
	if (nu > 100):
		return -99999999.9-nu
	if ((p < 0) | (p > 1)):
		return -99999999.9+p

	p = max(p, sys.float_info.epsilon)
	n_indv = len(events.keys())
	n_chr = len(events[0])
	assert n_chr == len(L), 'events does not match L in size'
	n_tot = n_indv * n_chr
	two_p = 2.0*p
	lk = 0.0

	log_1_minus_Gstar = np.zeros((n_indv, n_chr))
	log_1_minus_Gstar_escape = np.zeros((n_indv, n_chr))
	max_map = np.zeros(n_tot)
	all_events = [events[x][y] for x in events.keys() for y in events[x].keys()]

	tmp1 = np.zeros(n_chr)
	tmp2 = np.zeros(n_chr)
	for c in range(n_chr):
		tmp1[c] = np.log(1.0 - Gstar2(L[c], nu, p))
		tmp2[c] = np.log(1.0 - Gstar1(L[c], two_p))

	for c in range(n_chr):
		log_1_minus_Gstar[:,c] = tmp1[c]
		log_1_minus_Gstar_escape[:,c] = tmp2[c]

	n_partitions = np.array([np.power(3, len(x)) for x in all_events])
	uniq_partitions = np.unique(n_partitions)
	n_partitions_reshape = np.reshape(n_partitions, (n_indv, n_chr))

	idx = n_partitions_reshape == 1
	lk = lk + np.sum(2.0*log_1_minus_Gstar[idx] + log_1_minus_Gstar_escape[idx])

	for u in range(1, len(uniq_partitions)): # Could parallelize this loop
		n_part = uniq_partitions[u]
		width = len(dec2base3(n_part-1))
		partitions = np.array([dec2base3(x, width=width) for x in range(n_part)])
		idx_shape = np.where(n_partitions_reshape == n_part)
		idx = np.where(n_partitions == n_part)
		all_pos = np.array([events[idx_shape[0][x]][idx_shape[1][x]] for x in range(len(idx_shape[0]))])
		lka = np.tile(log_1_minus_Gstar[idx_shape], (n_part, 1))
		lkc = np.tile(log_1_minus_Gstar_escape[idx_shape], (n_part, 1))
		lkb = np.tile(np.zeros(len(idx[0])), (n_part, 1)) 
		max_map_idx = L[idx_shape[1]]
		I0 = partitions == 0
		I1 = partitions == 1
		I2 = partitions == 2
		
		# Only consider unique partitions
		uniq_I0 = np.vstack({tuple(row) for row in I0}).T
		uniq_I1 = np.vstack({tuple(row) for row in I1}).T
		uniq_I2 = np.vstack({tuple(row) for row in I2}).T
		
		uniq_I0 = uniq_I0[:,np.lexsort(np.rot90(uniq_I0.T))]
		uniq_I1 = uniq_I1[:,np.lexsort(np.rot90(uniq_I1.T))]
		uniq_I2 = uniq_I2[:,np.lexsort(np.rot90(uniq_I2.T))]

		lka_prime = lka[range(uniq_I0.shape[1]),:] 
		lkc_prime = lkc[range(uniq_I2.shape[1]),:] 
		sum_sites0 = np.sum(uniq_I0, axis=0)
		n_sites0 = uniq_I0.shape[1]
		sum_sites2 = np.sum(uniq_I2, axis=0)
		for n in range(1,n_sites0):
			J0 = sum_sites0 == n
			if (np.sum(J0) > 0):
				uniq_I00 = uniq_I0[:,J0]
				# Append boundaries
				uniq_I000 = np.concatenate((np.ones(uniq_I00.shape[1], dtype=bool)[:,np.newaxis], uniq_I00.T, np.ones(uniq_I00.shape[1], dtype=bool)[:,np.newaxis]), axis=1)
				all_pos0 = np.tile(np.concatenate((np.zeros(all_pos.shape[0])[:,np.newaxis], all_pos, max_map_idx[:,np.newaxis]), axis=1), (uniq_I000.shape[0], 1)).T
			
				uniq_I0000 =  np.reshape(np.tile(uniq_I000.T, (1, all_pos.shape[0])), all_pos0.shape)
				all_pos00 = np.ma.MaskedArray(all_pos0, ~uniq_I0000)
				x0 = np.diff([x.compressed() for x in all_pos00.T], axis=1)
				y0 = np.log(gstar2(x0[:,0], nu, p) * (1.0 - Fstar2(x0[:,-1], nu, p)))
				lka_prime[J0,:] = np.reshape(y0, lka_prime[J0,:].shape)
				if (n > 1):
					y1 = np.log(np.prod(fstar2(x0[:,1:-1], nu, p), axis=1))[:, np.newaxis]
					lka_prime[J0,:] += np.reshape(y1, lka_prime[J0,:].shape)

			J2 = sum_sites2 == n
			if (np.sum(J2) > 0):
				uniq_I22 = uniq_I2[:,J2]
				# Append boundaries
				uniq_I222 = np.concatenate((np.ones(uniq_I22.shape[1], dtype=bool)[:,np.newaxis], uniq_I22.T, np.ones(uniq_I22.shape[1], dtype=bool)[:,np.newaxis]), axis=1)
				all_pos2 = np.tile(np.concatenate((np.zeros(all_pos.shape[0])[:,np.newaxis], all_pos, max_map_idx[:,np.newaxis]), axis=1), (uniq_I222.shape[0], 1)).T

				uniq_I2222 = np.reshape(np.tile(uniq_I222.T, (1, all_pos.shape[0])), all_pos2.shape)
				all_pos22 = np.ma.MaskedArray(all_pos2, ~uniq_I2222)
				x2 = np.diff([x.compressed() for x in all_pos22.T], axis=1)
				y2 = np.log(gstar1(x2[:,0], two_p) * (1.0 - Fstar1(x2[:,-1], two_p)))
				lkc_prime[J2,:] = np.reshape(y2, lkc_prime[J2,:].shape)
				if (n > 1):
					y3 = np.log(np.prod(fstar1(x2[:,1:-1], two_p), axis=1))[:, np.newaxis]
					lkc_prime[J2,:] += np.reshape(y3, lkc_prime[J2,:].shape)
				
		# Assign prime matrices back to main matrix
		for counter, row in enumerate(uniq_I0.T):
			idx = np.where((I0 == row).all(axis=1))
			lka[idx,:] = lka_prime[counter, :]
			# No need to calculate lkb, as we can derive from lka by rearranging
			idx = np.where((I1 == row).all(axis=1))
			lkb[idx,:] = lka_prime[counter, :]
		for counter, row in enumerate(uniq_I2.T):
			idx = np.where((I2 == row).all(axis=1))
			lkc[idx,:] = lkc_prime[counter, :]

		tmp = lka + lkb + lkc
		ma = np.max(tmp, axis=0)
		tmp1 = ma + np.log(np.sum(np.exp(tmp - ma), axis=0))
		lk = lk + np.sum(tmp1)
	lk = lk + np.log(1.0 - p)
	print ("nu = %s p = %s lk = %s" % (nu, p, lk))
	return lk
