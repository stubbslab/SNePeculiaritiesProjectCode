import cantrips as can
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SNeFieldFits/'
    min_sn = 14
    grid_density = 24
    comoving_bin = 100

    rand_field_ids = [i  for i in range(1, 11)]
    #rand_field_ids = [i  for i in range(1, 101)]
    rand_files = ['OverdensityFit' + str(i) + '_MinNSN_' + str(min_sn) + '_GridDens_' + str(grid_density) + '_BinSize_' + str(comoving_bin) + '_Rand1_fits.txt'  for i in rand_field_ids]
    real_file = 'OverdensityFit' + str(1) + '_MinNSN_' + str(min_sn) + '_GridDens_' + str(grid_density) + '_BinSize_' + str(comoving_bin) + '_Rand0_fits.txt'

    best_fit_rand_vals = [np.min(np.array(can.readInColumnsToList(rand_file, file_dir = data_dir, n_ignore = 4, convert_to_float = 1, delimiter = ', ')[5]) / np.array(can.readInColumnsToList(rand_file, file_dir = data_dir, n_ignore = 4, convert_to_float = 1, delimiter = ', ')[6])) for rand_file in rand_files]
    sorted_fit_rand_vals, sorted_rand_field_ids = can.safeSortOneListByAnother(best_fit_rand_vals, [best_fit_rand_vals, rand_field_ids])
    print ('Here are the randomized field ids, in order of best (most improved) fit: ')
    print(sorted_rand_field_ids)
    best_fit_real_val = np.min(np.array(can.readInColumnsToList(real_file, file_dir = data_dir, n_ignore = 4, convert_to_float = 1, delimiter = ', ')[5]) / np.array(can.readInColumnsToList(real_file, file_dir = data_dir, n_ignore = 4, convert_to_float = 1, delimiter = ', ')[6]))
    n_better = len([val for val in best_fit_rand_vals if val > best_fit_real_val])

    n_hist_bins = 20
    hist = plt.hist(best_fit_rand_vals, bins = n_hist_bins, color = 'white', ec = 'k')
    hist_levels = hist[0]
    hist_centers= hist[1]
    hist_max = np.max(hist_levels)
    max_center = np.max(hist_centers)
    plt.axvline(best_fit_real_val, c = 'r', linestyle = '--')
    plt.text(best_fit_real_val, hist_max / 2, 'Real improvement over null', rotation = 90, c = 'r', horizontalalignment = 'right', verticalalignment = 'center')
    plt.text(max_center, hist_max * 9 / 10, 'Real data < ' + str(n_better) + '/' + str(len(best_fit_rand_vals)) + ' of bootstrapped data', c = 'k', horizontalalignment = 'right', verticalalignment = 'center')
    plt.xlabel('BEST improvement over null fit')
    plt.ylabel(r'$N$ permutation tests')
    plt.title('Permutation test of SNe correlations')
    plt.savefig('/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/' + 'HistOfNullImprovements_GridDensity' + str(grid_density) + '_MinSN' + str(min_sn) + '_ComovingBin' + str(comoving_bin) + '.pdf')
    plt.show()
