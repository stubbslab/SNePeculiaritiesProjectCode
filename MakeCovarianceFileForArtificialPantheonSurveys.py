import cantrips as can

"""
Take in the Pantheon+ SNe covariance file, and add new lines for new
   artificial SNe.  These new lines are just 0's (assume no
   covariances with artificial SNe).  The 0's just need to be placed
   at the right spots to preserve the actual covariances between the
   real data points.
"""

n_extra_sn = 300

base_covariance_file = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/OriginalSNDataFiles/Pantheon+SH0ES_122221_1.cov'
covariance_file_to_save = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_covariances_extra' + str(n_extra_sn) + '.txt'

orig_covariance = [float(elem) for elem in can.readInColumnsToList(base_covariance_file, verbose = 0)[0]]
orig_n_sn = int(orig_covariance[0])
orig_covariance = orig_covariance[1:]
orig_lines = [orig_covariance[i * orig_n_sn:(i+1) * orig_n_sn] for i in range(orig_n_sn)]
new_lines = [orig_line + [0.0 for i in range(n_extra_sn)] for orig_line in orig_lines]
new_lines = new_lines + [[0.0 for i in range(n_extra_sn + orig_n_sn)] for j in range(n_extra_sn)]

new_covariances = can.flattenListOfLists(new_lines)
can.saveListsToColumns([int(orig_n_sn + n_extra_sn)] + new_covariances, covariance_file_to_save, '')
