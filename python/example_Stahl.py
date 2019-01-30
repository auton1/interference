import simStahl as sim
import stahlLogLk as lk
import numpy as np
import scipy.optimize as op
import pdb

# Run a simulation study of the Stahl model (non-quad / phase known data).
#
# Usage: [nu_est, p_est] = example_Stahl(nu, p, cM_map_len, N_indv)
#
# Output:
# 	nu_est: estimate of interference parameter
#	p_est: estimate of escape parameter
#	lk_max: log likelihood of fitted model
#   covariance_matrix: covariance matrix estimated as inverse of fisher information matrix.
#
# Input:
#	nu: interference parameter to simulate
#   p: escape parameter to simulate
#	cM_map_len:	vector of chromosome map lengths in centiMorgans
#	N_indv: number of individuals to simulate
#
def opt_fun(x, args):
	return -1.0 * lk.stahlLogLk(args[0], args[1], x[0], x[1])

def opt_fun_quad(x, args):
	return -1.0 * lk.stahlLogLk_quad(args[0], args[1], x[0], x[1])

def example_Stahl(nu, p, cM_map_len = np.array([200, 250]), n_indv = 300):
	assert ((nu >= 0.1) & (nu <= 50.0)), 'Require nu in range 0.1 to 50'
	assert ((p >= 0.0) & (p <= 1.0)), 'Require p in range 0 to 1'
	
	L = cM_map_len / 100.0 # Convert to Morgans
	events = sim.simStahl(n_indv, L, nu, p)
	
	lk_truth = lk.stahlLogLk(events, L, nu, p)
	print("lk of truth = %s" % (lk_truth)) 
	
	x0 = [0.1, 0.1]
	#x0 = [nu, p]
	res = op.minimize(opt_fun, x0, args=[events, L], method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
	nu_est = res.x[0]
	p_est = res.x[1]
	return([nu_est, p_est])

def example_Stahl_quad(nu, p, cM_map_len = np.array([200, 250]), n_indv = 300):
	assert ((nu >= 0.1) & (nu <= 50.0)), 'Require nu in range 0.1 to 50'
	assert ((p >= 0.0) & (p <= 1.0)), 'Require p in range 0 to 1'
	
	L = cM_map_len / 100.0 # Convert to Morgans
	events = sim.simStahl_quad(n_indv, L, nu, p)
	
	lk_truth = lk.stahlLogLk_quad(events, L, nu, p)
	print("lk of truth = %s" % (lk_truth)) 
	
	x0 = [0.1, 0.1]
	#x0 = [nu, p]
	res = op.minimize(opt_fun_quad, x0, args=[events, L], method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
	nu_est = res.x[0]
	p_est = res.x[1]
	return([nu_est, p_est])
	
if __name__ == "__main__":
	[nu_est, p_est] = example_Stahl(7.0, 0.05)
	[nu_est_quad, p_est_quad] = example_Stahl_quad(7.0, 0.05)
