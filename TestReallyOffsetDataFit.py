import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import scipy

def log_probability_no_offset(theta, xs0, ys0, xs1, ys1):
    #print ('theta = ' + str(theta))
    log_prob = -0.5 * (np.sum((np.array(ys0) - single_funct(np.array(xs0), theta[0], theta[2], theta[3])) ** 2.0) + np.sum((np.array(ys1) - single_funct(np.array(xs1), theta[1], theta[2], theta[3])) ** 2.0) )
    #print ('log_prob = ' + str(log_prob))

    return log_prob

if __name__ == "__main__":

     n_samples = 100
     xs0 = [np.random.random() * 10 for i in range(n_samples // 2)]
     xs1 = [np.random.random() * 10 for i in range(n_samples // 2)]
     a0, a1, b, c = [0, 10, -2, 0.5]
     single_funct = lambda xs, a, b, c: a + b * xs + c * xs ** 2.0
     ys0 = [single_funct(elem, a0, b, c) + np.random.random() * 2 for elem in xs0]
     ys1 = [single_funct(elem, a1, b, c) + np.random.random() * 2 for elem in xs1]

     init_no_offset_fit = scipy.optimize.curve_fit( lambda xs, a, b, c: single_funct(np.array(xs), a, b, c), xs0 + xs1, ys0 + ys1 )

     pos_no_offset = init_no_offset_fit[0] + 1e-4 * np.random.randn(32, len(init_no_offset_fit[0]) )
     nwalkers, ndim = pos_no_offset.shape

     loaded_log_probability_no_offset = lambda theta: log_probability_no_offset([theta[0], theta[0], theta[1], theta[2]], xs0, ys0, xs1, ys1)
     sampler_no_offset = emcee.EnsembleSampler(nwalkers, ndim, loaded_log_probability_no_offset, args=())
     sampler_no_offset.run_mcmc(pos_no_offset, 5000, progress = True )
     print ('np.shape(pos_no_offset) = ' + str(np.shape(pos_no_offset)))
     flat_samples_no_offset = sampler_no_offset.get_chain(discard=100, thin=15, flat=True)


     init_with_offset_fit = scipy.optimize.curve_fit( lambda xs, a0, b, c, a1: single_funct(np.array(xs[0:len(xs) // 2]), a0, b, c).tolist() + single_funct(np.array(xs[len(xs) // 2:]), a1, b, c).tolist(), xs0 + xs1, ys0 + ys1 )
     pos_with_offset = init_with_offset_fit[0] + 1e-4 * np.random.randn(32, len(init_with_offset_fit[0]) )
     nwalkers, ndim = pos_with_offset.shape
     loaded_log_probability_with_offset = lambda theta: log_probability_no_offset([theta[0], theta[3], theta[1], theta[2]], xs0, ys0, xs1, ys1)
     print ('[nwalkers, ndim ] = ' + str([nwalkers, ndim]))
     print ('np.shape(pos_with_offset) = ' + str(np.shape(pos_with_offset)))
     sampler_with_offset = emcee.EnsembleSampler(nwalkers, ndim , loaded_log_probability_with_offset, args=())
     sampler_with_offset.run_mcmc(pos_with_offset, 5000, progress = True )
     flat_samples_with_offset = sampler_with_offset.get_chain(discard=100, thin=15, flat=True)


     plt.scatter(xs0, ys0, c = 'k')
     plt.scatter(xs1, ys1, c = 'g')
     plt.show()

     guess_truths = [(a0 + a1)/2, b, c]
     print ('guess_truths = ' + str(guess_truths))
     mcmc_truths = [ np.median([flat_samples_no_offset[i][j] for i in range(len(flat_samples_no_offset))]) for j in range(len(flat_samples_no_offset[0])) ]
     mcmc_no_offset_stds = [ np.std([flat_samples_no_offset[i][j] for i in range(len(flat_samples_no_offset))]) for j in range(len(flat_samples_no_offset[0])) ]
     mcmc_with_offset_stds = [ np.std([flat_samples_with_offset[i][j] for i in range(len(flat_samples_with_offset))]) for j in range(len(flat_samples_with_offset[0])) ]
     print ('mcmc_truths = ' + str(mcmc_truths))
     print ('init_no_offset_fit[0] = ' + str(init_no_offset_fit[0]))
     fig_no_offset = corner.corner(flat_samples_no_offset, labels=['a0', 'b','c'], truths=guess_truths)
     plt.show()
     fig_with_offset = corner.corner(flat_samples_with_offset, labels=['a0','b','c',  'a1'], truths=guess_truths + [guess_truths[0]])
     plt.show()

     print ('mcmc_no_offset_stds = ' + str(mcmc_no_offset_stds ))
     print ('mcmc_with_offset_stds = ' + str(mcmc_with_offset_stds ))
