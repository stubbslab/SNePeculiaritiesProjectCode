import loadSN as lsn
import cantrips as c
import numpy as np
import CosmologicalParameterArchive as cpa
import binData as bd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib as mpl

def plotFinalResidualsInTimeSpaceOfRingermacherResultsOld(rand_ringermacher_replicators, true_ringermacher_replicator = None, n_fig_units = [7,5], figsize = [8,12],
                                                       true_ringermacher_indeces = [4,3], x_lims = [0.39, 1.01], y_lims = [-0.24, 0.24], rand_color = 'k', true_color = 'k',
                                                       save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/RingermacherResponse/Plots/', save_fig_name = 'SampleRandomPlots_True_at_4_3.pdf',):
    f, axarr = plt.subplots(*n_fig_units, figsize = figsize, squeeze = False, sharex = True, sharey = True)
    for i in range(n_fig_units[0]):
        if i == n_fig_units[0] - 1:
            xlabel_to_plot = xlabel
        else:
            xlabel_to_plot = ''
        for j in range(n_fig_units[1]):
            if j == 0:
                ylabel_to_plot = ylabel
            else:
                ylabel_to_plot = ''
            if not(true_ringermacher_replicator is None) and [i,j] == true_ringermacher_indeces:
                color = true_color
                rr_to_plot = true_ringermacher_replicator
            else:
                color = rand_color
                rr_to_plot = rand_ringermacher_replicators[i * n_fig_units[1] + j]
            time_cut_indeces = rr_to_plot.time_cut_indeces
            #rr_to_plot.addSinglePlotToArray(axarr[i,j], rr_to_plot.tauMeas_bins[time_cut_indeces[0]:time_cut_indeces[1]], rr_to_plot.smooth_aDerivs_diff[time_cut_indeces[0]:time_cut_indeces[1]], x_lims, y_lims, legend, xlabel_to_plot, ylabel_to_plot, '', color = color, labelsize = labelsize, xticks = xticks, yticks = yticks)
            rr_to_plot.addSinglePlotToArray(axarr[i,j], rr_to_plot.tauMeas_bins[time_cut_indeces[0]:time_cut_indeces[1]], rr_to_plot.smooth_aDerivs_diff[time_cut_indeces[0]:time_cut_indeces[1]], x_lims, y_lims, legend, xlabel_to_plot, ylabel_to_plot, '', color = color, labelsize = labelsize, yticks = yticks)
    mpl.rcParams['xtick.labelsize'] = tick_label_size
    mpl.rcParams['ytick.labelsize'] = tick_label_size
    plt.subplots_adjust(hspace = 0.1, wspace = 0.15)
    axarr[0,n_fig_units[1] // 2].set_title(title, fontsize = titlesize)
    plt.tight_layout()
    #f.suptitle(title)
    if save_fig:
        plt.savefig(save_dir + save_fig_name)
    if show_fig:
       plt.show()
    plt.close()


def plotFinalResidualsInTimeSpaceOfRingermacherResults(rand_ringermacher_replicators, true_ringermacher_replicator = None, n_fig_units = [7,5], figsize = [8,12],
                                                       true_ringermacher_indeces = [4,2], x_lims = [[0.0, 2.3], [0.39, 1.01]], y_lims = [[-0.7, 0.7], [-0.24, 0.24]], rand_colors = ['k','k'], true_color = 'k',
                                                       ylabels = [r'$\Delta \mu$ [mag]' , r'$\Delta G(d{\overline{\Delta a}}/dt)$'], xlabels = ['z',r'$t$'],
                                                       titles = ['Pantheon-like Data', 'R20 Analyisis of Pantheon-like Data'], suptitle = r'R20 Processing Applied to Artificial Pantheon-like $\Lambda$CDM Data', titlesize = 20,
                                                       legend = None, labelsize = 12, labelpad = -5, tick_label_size = 10.0,
                                                       save_fig = 0, show_fig = 1, xticks = [0.5, 0.7, 0.9], yticks = [-0.2, -0.1, 0.0, 0.1, 0.2],
                                                       save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/RingermacherResponse/Plots/', save_fig_name = 'SampleRandomPlots_True_at_4_3.pdf',):
    f, axarr = plt.subplots(n_fig_units[0], n_fig_units[1] * 2, figsize = figsize, squeeze = False, sharex = False, sharey = False)

    for i in range(n_fig_units[0]):
        if i == n_fig_units[0] - 1:
            xlabel_to_plot = xlabels[1]
            title_to_plot = ''
        elif i == 0:
            title_to_plot = titles[1]
            xlabel_to_plot = ''
        else:
            title_to_plot = ''
            xlabel_to_plot = ''
        #Plot the processed data in one column...
        for j in range(n_fig_units[1]):
            if not(true_ringermacher_replicator is None) and [i,j] == true_ringermacher_indeces:
                color = true_color
                rr_to_plot = true_ringermacher_replicator
            else:
                color = rand_colors[1]
                rr_to_plot = rand_ringermacher_replicators[i * n_fig_units[1] + j]
            time_cut_indeces = rr_to_plot.time_cut_indeces
            #rr_to_plot.addSinglePlotToArray(axarr[i,j], rr_to_plot.tauMeas_bins[time_cut_indeces[0]:time_cut_indeces[1]], rr_to_plot.smooth_aDerivs_diff[time_cut_indeces[0]:time_cut_indeces[1]], x_lims, y_lims, legend, xlabel_to_plot, ylabel_to_plot, '', color = color, labelsize = labelsize, xticks = xticks, yticks = yticks)
            rr_to_plot.addSinglePlotToArray(axarr[i,j*2+1], rr_to_plot.tauMeas_bins[time_cut_indeces[0]:time_cut_indeces[1]], rr_to_plot.smooth_aDerivs_diff[time_cut_indeces[0]:time_cut_indeces[1]], x_lims[1], y_lims[1], legend, xlabel_to_plot, ylabels[1], title_to_plot, color = color, labelsize = labelsize, yticks = yticks, labelpad = labelpad)
            if i != n_fig_units[0] - 1:
                axarr[i,j*2+1].set_xticklabels('')
            axarr[i,j*2+1].tick_params(axis = 'x', direction = 'in', pad = 2)
            axarr[i,j*2+1].tick_params(axis = 'y', direction = 'in', pad = 2)
        #... and plot the original z vs \Delta \mu data in the other .
        if i == n_fig_units[0] - 1:
            axarr[i, j*2].set_xlabel(xlabels[0], labelpad = labelpad)
        if i == 0:
            axarr[i,j*2].set_title(titles[0])
        for j in range(n_fig_units[1]):
            if not(true_ringermacher_replicator is None) and [i,j] == true_ringermacher_indeces:
                color = true_color
                rr_to_plot = true_ringermacher_replicator
            else:
                color = rand_colors[0]
                rr_to_plot = rand_ringermacher_replicators[i * n_fig_units[1] + j]
            axarr[i,j*2].scatter(rr_to_plot.sorted_zs, rr_to_plot.sorted_muDiffs, marker = '.', s = 5, c = color)
            #axarr[i,j*2].errorbar(rr_to_plot.sorted_zs, rr_to_plot.sorted_muDiffs, yerr = rr_to_plot.sorted_muErrs, fmt = 'none', color = color)
            axarr[i,j*2].set_xlim(x_lims[0])
            axarr[i,j*2].set_ylim(y_lims[0])
            #rr_to_plot.addSinglePlotToArray(axarr[i,j], rr_to_plot.tauMeas_bins[time_cut_indeces[0]:time_cut_indeces[1]], rr_to_plot.smooth_aDerivs_diff[time_cut_indeces[0]:time_cut_indeces[1]], x_lims, y_lims, legend, xlabel_to_plot, ylabel_to_plot, '', color = color, labelsize = labelsize, xticks = xticks, yticks = yticks)
            #rr_to_plot.addSinglePlotToArray(axarr[i,j], rr_to_plot.tauMeas_bins[time_cut_indeces[0]:time_cut_indeces[1]], rr_to_plot.smooth_aDerivs_diff[time_cut_indeces[0]:time_cut_indeces[1]], x_lims, y_lims, legend, xlabel_to_plot, ylabel_to_plot, '', color = color, labelsize = labelsize, yticks = yticks)
            axarr[i, j*2].set_ylabel(ylabels[0], labelpad = labelpad)
            if i != n_fig_units[0] - 1:
                axarr[i,j*2].set_xticklabels('')
            axarr[i,j*2].tick_params(axis = 'x', direction = 'in', pad = 2)
            axarr[i,j*2].tick_params(axis = 'y', direction = 'in', pad = 2)
    mpl.rcParams['xtick.labelsize'] = tick_label_size
    mpl.rcParams['ytick.labelsize'] = tick_label_size

    plt.subplots_adjust(hspace = 0.1, wspace = 0.15)
    #axarr[0,n_fig_units[1] // 2].set_title(title, fontsize = titlesize)
    #plt.tight_layout()
    f.suptitle(suptitle, fontsize = titlesize, y = 0.92)
    if save_fig:
        plt.savefig(save_dir + save_fig_name)
    if show_fig:
       plt.show()
    plt.close()

def generateHistOfChosenFrequency(rand_ringermacher_replicators, chosen_frequency, true_ringermacher_replicator = None, bins = 21,
                                  save_fig = 0, show_fig = 1, save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/RingermacherResponse/Plots/', save_fig_name = 'RingermacherHistogramAtChosenFrequency.pdf',
                                  xlabel = r'Power in $\Delta G(d{\overline{\Delta a}}/dt)$ at ', ylabel = r'Number of randomizations', title = 'Power of Artificial Pantheon Data at ',
                                  show_true_rr = 1, textsize = 12, labelsize = 14, titlesize = 16, color = 'k', figsize = [10, 8]):

    if true_ringermacher_replicator is None:
        frequencies = rand_ringermacher_replicators[0].frequencies
    else:
        frequencies = true_ringermacher_replicator.frequencies
    freq_index = np.argmin(np.abs(np.array(frequencies) - chosen_frequency))
    rand_power_at_freq = [rr.fft_power[freq_index] for rr in rand_ringermacher_replicators]
    f = plt.figure(figsize=figsize)
    (n, bins, patches) = plt.hist(rand_power_at_freq, bins = bins, color = 'white', edgecolor = color,)

    if not(true_ringermacher_replicator is None):
        true_power_at_freq = true_ringermacher_replicator.fft_power[freq_index]
    plt.axvline(true_power_at_freq, c = 'k')
    plt.xlabel(xlabel + str(chosen_frequency) + r'$H$Hz', fontsize = labelsize)
    plt.ylabel(ylabel, fontsize = labelsize)
    plt.title(title + str(chosen_frequency) + r'$H$Hz', fontsize = titlesize)

    frac_of_rands_above_true_power = len([elem for elem in rand_power_at_freq if elem > true_power_at_freq]) / len(rand_power_at_freq)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #textstr = str(c.round_to_n(len([power for power in rand_power_at_freq if power > true_power_at_freq]) / len(rand_ringermacher_replicators) * 100,3)) + '\% of random realization of $\Lambda$CDM ' + ' Pantheon-like data ' + '\n have more power at ' + str(chosen_frequency) + r'$H$Hz' + ' than the \n real Pantheon data.'
    textstr = str(c.round_to_n(frac_of_rands_above_true_power * 100, 3)) + r'$\%$ of random realizations of $\Lambda$CDM ' + '\nPantheon-like data ' + ' have more power at \n' + str(chosen_frequency) + r'$H$Hz' + ' than the real Pantheon data.'
    #textstr = str(len([power for power in rand_power_at_freq if power > true_power_at_freq])) + '/' + str(len(rand_ringermacher_replicators))
    plt.text(bins[int(len(bins) * 0.55)],  max(n) * 0.9, textstr, fontsize=textsize,
             verticalalignment='top', bbox=props)

    plt.tight_layout()
    if save_fig:
        plt.savefig(save_dir + save_fig_name)
    if show_fig:
        plt.show()
    plt.close()


def showDistributionOfRingermacherPeriodogram(rand_ringermacher_replicators, true_ringermacher_replicator = None, save_fig = 0, show_fig = 1, save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/RingermacherResponse/Plots/', save_fig_name = 'RingermacherFourierPlot.pdf',
                                              xlabel = r'Normalized Cosmic time frequency, $H$Hz', ylabel = r'Power of Oscillation in $\Delta G(d{\overline{\Delta a}}/dt)$', title = 'Fourier Spectra of Real and Artificial Data',
                                              show_true_rr = 1, labelsize = 14, titlesize = 17, vert_line_freqs = [7.5]):

    n_rand = len(rand_ringermacher_replicators)
    fft_amps = [rr.fft_power for rr in rand_ringermacher_replicators]
    print ('len(rand_ringermacher_replicators[0].fft_amplitudes) = ' + str(len(rand_ringermacher_replicators[0].fft_power)))
    if true_ringermacher_replicator is None:
        frequencies = rand_ringermacher_replicators[0].frequencies
    else:
        frequencies = true_ringermacher_replicator.frequencies
    means = np.mean(fft_amps, axis = 0)
    stds = np.std(fft_amps, axis = 0)
    print ('len(means) = ' + str(len(means)))
    #plt.plot(means)
    #plt.plot(means + 1.0 * stds)
    #plt.plot(means - 1.0 * stds)
    percent_bounds = [0.0, 0.5, 0.8, 0.9, 0.95, 0.99]
    colors = ['g', 'r', 'b', 'pink', 'yellow']
    fft_amp_by_frequencies = np.transpose(fft_amps)
    stacked_fft_amp_by_frequencies = [sorted(fft_amp_by_freq) for fft_amp_by_freq in fft_amp_by_frequencies]
    bounds_plots = [plt.fill_between(frequencies, [stacked_fft_amp[int(len(stacked_fft_amp) * percent_bounds[i+1])] for stacked_fft_amp in stacked_fft_amp_by_frequencies], [stacked_fft_amp[int(len(stacked_fft_amp) * percent_bounds[i])] for stacked_fft_amp in stacked_fft_amp_by_frequencies], alpha = 0.3, color = colors[i]) for i in range(len(percent_bounds) - 1)]
    #oneSigBound = plt.fill_between(frequencies, [max(val, 0.0) for val in (means - 1.0 * stds)], means + 1.0 * stds, alpha = 0.3, color = 'r')
    #twoSigBound = plt.fill_between(frequencies, [max(val, 0.0) for val in (means - 2.0 * stds)], means + 2.0 * stds, alpha = 0.2, color = 'b')
    plt.xlabel(xlabel, fontsize = labelsize)
    plt.ylabel(ylabel, fontsize = labelsize)
    plt.title(title, fontsize = titlesize)
    plt.tight_layout()
    for vert_line_freq in vert_line_freqs:
        plt.axvline(vert_line_freq, color = 'grey', alpha = 0.5)
    if not (true_ringermacher_replicator is None) and show_true_rr:
        true_plot = plt.plot(frequencies, true_ringermacher_replicator.fft_power, c = 'k')[0]
        plt.legend([*bounds_plots, true_plot], [str(int(frac * 100)) + r"$\%$" + r' of randomizations have less Fourier power' for frac in percent_bounds[1:]] + ['True data'])
    else:
        plt.legend(['oneSigBound, twoSigBound'])
    if save_fig:
        plt.savefig(save_dir + save_fig_name)
    if show_fig:
       plt.show()
    plt.close()


class RingermacherReplicator:

    def addSinglePlotToArray(self, ax, x_data, y_data, xlims, ylims, legend_text, xlabel, ylabel, title, color = 'k', labelsize = 12, xticks = None, yticks =  None, labelpad = 0):
        scat = ax.scatter(x_data, y_data, marker = '.', c = color)
        #scat = ax.scatter(x_data[:], y_data[:], marker = '.')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        if not (legend_text is None):
            ax.legend([scat], [legend_text])
        ax.set_xlabel(xlabel, fontsize = labelsize, labelpad = labelpad)
        ax.set_ylabel(ylabel, fontsize = labelsize, labelpad = labelpad)
        if not(xticks is None):
            ax.set_yticks(xticks)
        if not(yticks is None):
            ax.set_yticks(yticks)
        ax.set_title(title)
        return 1

    def showFourierTransform(self, save_fig = 0, show_fig = 1, save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/RingermacherResponse/Plots/', save_fig_name = 'RingermacherFourierPlot.pdf',
                             xlabel = r'Normalized Cosmic time frequency, $H$Hz', ylabel = r'Amplitude of Oscillation in $da/dt$', title = 'Fourier Spectrum of Final Plot'):
        plt.rc('font', family='serif')
        plt.rc('text', usetex=True)
        f, axarr = plt.subplots(1,1,squeeze = False)
        axarr[0,0].plot(self.fft_power)
        axarr[0,0].set_xlabel(xlabel)
        axarr[0,0].set_ylabel(ylabel)
        axarr[0,0].set_title(title)
        if save_fig:
            plt.savefig(save_dir + save_fig_name)
        if show_fig:
           plt.show()
        plt.close()


    def showDataSteps(self, show_resids = 1, show_binned_a_resids = 1, show_derivs_of_binned_resids = 1, show_high_pass_smooth = 1, show_low_pass_smooth = 1, show_smooth_diff = 1,
                      shared_tau_lims = [0.4, 1.0], a_lims = [-0.07, 0.07], a_dot_lims = [-0.2, 0.2], a_dot_ticks = [-0.2, -0.1, 0.0, 0.1, 0.2], vertical_plot_elem_size = 1.3, horizontal_plot_size = 5, xlabel = r'Normalized Cosmic time, $t$',
                       show_fig = 1, save_fig = 0, save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/RingermacherResponse/Plots/', save_fig_name = 'RingermacherPlot.pdf', title = 'Processing Steps of R20',
                       titlesize = 18 ):
        plt.rc('font', family='serif')
        plt.rc('text', usetex=True)
        n_plots = sum([show_resids, show_binned_a_resids, show_derivs_of_binned_resids, show_high_pass_smooth, show_low_pass_smooth, show_smooth_diff, ])
        figsize = (horizontal_plot_size, vertical_plot_elem_size * n_plots)
        f, axarr = plt.subplots(n_plots, 1, figsize = figsize, sharex = True, squeeze = False)
        current_plot = 0
        if show_resids:
            y_lims = a_lims
            ylabel = r'$\Delta a$'
            self.addSinglePlotToArray(axarr[current_plot, 0], self.sorted_tauMeas, self.sorted_asMeasResiduals, shared_tau_lims, y_lims, r'$a$ of $t$ residuals', '', ylabel, '')
            current_plot = current_plot + 1
        if show_binned_a_resids:
            y_lims = a_lims
            ylabel = r'$\overline{\Delta a}$'
            self.addSinglePlotToArray(axarr[current_plot, 0], self.tauMeas_bins[self.time_cut_indeces[0]:self.time_cut_indeces[1]], self.binned_aMeasResiduals[self.time_cut_indeces[0]:self.time_cut_indeces[1]], shared_tau_lims, y_lims, r'binned $a$ of $t$ residuals', '', ylabel, '')
            current_plot = current_plot + 1
        if show_derivs_of_binned_resids:
            y_lims = a_dot_lims
            ylabel = r'$d{\overline{\Delta a}}/dt$'
            self.addSinglePlotToArray(axarr[current_plot, 0], self.tauMeas_bins[self.time_cut_indeces[0]:self.time_cut_indeces[1]], self.binned_aMeasDerivs[self.time_cut_indeces[0]:self.time_cut_indeces[1]], shared_tau_lims, y_lims, r'derivs of binned $a$ of $t$ residuals', '', ylabel, '', yticks = a_dot_ticks)
            current_plot = current_plot + 1
        if show_high_pass_smooth:
            y_lims = a_dot_lims
            ylabel = r'$G_{0.05}(d{\overline{\Delta a}}/dt)$'
            self.addSinglePlotToArray(axarr[current_plot, 0], self.tauMeas_bins[self.time_cut_indeces[0]:self.time_cut_indeces[1]], self.high_pass_smooth_aDerivs[self.time_cut_indeces[0]:self.time_cut_indeces[1]], shared_tau_lims, y_lims, r'high-pass smoothed derivs', '', ylabel, '', yticks = a_dot_ticks)
            current_plot = current_plot + 1
        if show_low_pass_smooth:
            y_lims = a_dot_lims
            ylabel = r'$G_{0.13}(d{\overline{\Delta a}}/dt)$'
            self.addSinglePlotToArray(axarr[current_plot, 0], self.tauMeas_bins[self.time_cut_indeces[0]:self.time_cut_indeces[1]], self.low_pass_smooth_aDerivs[self.time_cut_indeces[0]:self.time_cut_indeces[1]], shared_tau_lims, y_lims, r'low-pass smoothed derivs', '', ylabel, '', yticks = a_dot_ticks)
            current_plot = current_plot + 1
        #if show_smooth_diff:
        #    y_lims = a_dot_lims
        #    ylabel = r'$\dot{\overline{\Delta a}}^{\Delta G_{k}}$'
        #    print ('np.mean( self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] ) = '+ str(np.mean( self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] )))
        #    self.addSinglePlotToArray(axarr[current_plot, 0], self.tauMeas_bins[self.time_cut_indeces[0]:self.time_cut_indeces[1]], self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] - np.mean( self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] ), shared_tau_lims, y_lims, r'high-pass $-$ low-pass', '', ylabel, '', color = 'r')
        #    #current_plot = current_plot + 1
        if show_smooth_diff:
            y_lims = a_dot_lims
            ylabel = r'$\Delta G(d{\overline{\Delta a}}/dt)$'
            print ('np.mean( self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] ) = '+ str(np.mean( self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] )))
            self.addSinglePlotToArray(axarr[current_plot, 0], self.tauMeas_bins[self.time_cut_indeces[0]:self.time_cut_indeces[1]], self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]] , shared_tau_lims, y_lims, r'high-pass $-$ low-pass', '', ylabel, '', color = 'k', yticks = a_dot_ticks)
            current_plot = current_plot + 1
        axarr[-1,0].set_xlabel(xlabel)
        axarr[0,0].set_title(title, fontsize = titlesize)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0.2)
        if save_fig:
            plt.savefig(save_dir + save_fig_name)
        if show_fig:
           plt.show()
        plt.close()

    def quantifyOscillationStrength(self):
        mean_aDeriv = np.mean(self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]])
        mean_sub_aDerivs = [0.0 if np.isnan(aDeriv) else aDeriv - mean_aDeriv for aDeriv in self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]]  ]
        n_bins = len(mean_sub_aDerivs)
        fourier_amplitude_from_python = np.abs(np.fft.rfft(mean_sub_aDerivs))
        my_cos_sum = np.array([np.sum([mean_sub_aDerivs[i] * np.cos(2.0 * np.pi / n_bins * k * i) for i in range(len(mean_sub_aDerivs))]) for k in range(0, n_bins // 2)])
        my_sin_sum = np.array([np.sum([mean_sub_aDerivs[i] * np.sin(2.0 * np.pi / n_bins * k * i) for i in range(len(mean_sub_aDerivs))]) for k in range(0, n_bins // 2)])
        fourier_power_from_me = (my_cos_sum ** 2.0 + my_sin_sum ** 2.0 ) / self.n_bins
        #self.fft_of_binned_adot = np.fft.rfft([0.0 if np.isnan(aDeriv) else aDeriv - mean_aDeriv for aDeriv in self.smooth_aDerivs_diff[self.time_cut_indeces[0]:self.time_cut_indeces[1]]  ])

        self.python_fft_amplitudes = fourier_amplitude_from_python
        self.my_fft_power = fourier_power_from_me
        self.fft_power = self.my_fft_power
        time_range = self.time_cuts[1] - self.time_cuts[0]
        self.frequencies = [1.0 / time_range * i for i in range(len(self.fft_power))]

        #plt.plot(np.sqrt(np.real(self.fft_of_binned_adot) ** 2.0 + np.imag(self.fft_of_binned_adot) ** 2.0))
        #plt.show()

    def __init__(self, OmM = 0.27, OmL = 0.73, H0 = 68.75, Omk = 0.0, rand_around_data = 0, rand_around_canon = 0, time_cuts = [0.46, 1.0], n_bins = 128, n_t_interps = 1001):
        self.time_cuts = time_cuts
        self.OmM, self.OmL, self.H0, self.Omk = [OmM, OmL, H0, Omk]
        all_sn = lsn.loadSN(1, pull_extinctions = 0, OmL = OmL, OmM = OmM, H0 = H0, zHD = 1)
        print ('len(all_sn) = ' + str(len(all_sn)))
        all_zs = [sn['z'] for sn in all_sn]
        all_muMeas = [sn['mu'] for sn in all_sn]
        all_muErrs = [sn['muErr'] for sn in all_sn]
        all_muDiffs= [sn['muDiff'] for sn in all_sn]
        all_muTheory = np.array(all_muMeas) - np.array(all_muDiffs)
        if rand_around_data: all_muMeas = [np.random.normal(all_muMeas[i], all_muErrs[i]) for i in range(len(all_muMeas))]
        if rand_around_canon: all_muMeas = [np.random.normal(all_muTheory[i], all_muErrs[i]) for i in range(len(all_muTheory))]
        all_muDiffs = np.array(all_muMeas) - all_muTheory
        all_dLs = [sn['dl'] for sn in all_sn]
        self.cosmo_arch = cpa.CosmologicalParameterArchive()
        self.speed_of_light = self.cosmo_arch.getc()

        self.dH = self.speed_of_light / self.H0


        self.sorted_zs, self.sorted_muMeas, self.sorted_muErrs, self.sorted_muDiffs, self.sorted_muTheory, self.sorted_dLs = c.safeSortOneListByAnother(all_zs, [all_zs, all_muMeas, all_muErrs, all_muDiffs, all_muTheory, all_dLs])

        dl_of_mu = lambda mus: 10.0 ** ((np.array(mus) - 25.0)/ 5.0)
        self.sorted_dLTheory = dl_of_mu(self.sorted_muTheory)
        self.sorted_dLMeas = dl_of_mu(self.sorted_muMeas)

        a_of_z = lambda zs: 1.0 / (1.0 + np.array(zs))
        self.sorted_as = a_of_z(self.sorted_zs)

        Ys_of_funct = lambda a_s, dLs, dH:  np.array(a_s) * np.array(dLs) / dH
        self.sorted_YTheory = Ys_of_funct(self.sorted_as, self.sorted_dLTheory, self.dH)
        self.sorted_YMeas = Ys_of_funct(self.sorted_as, self.sorted_dLMeas, self.dH)

        aDeltaYs_of_funct = lambda a_s, Ys: np.array([a_s[i] * (Ys[i] - Ys[i-1]) for i in range(1, len(a_s))])

        self.sorted_aDeltaYsTheory = aDeltaYs_of_funct(self.sorted_as, self.sorted_YTheory)
        self.sorted_aDeltaYsMeas = aDeltaYs_of_funct(self.sorted_as, self.sorted_YMeas)

        taus_of_funct = lambda a_s, Ys, sub_val: np.array([1.0 - np.sum([a_s[i] * (Ys[i] - Ys[i-1]) for i in range(1, j+1)]) for j in range(1, len(a_s))]) - sub_val
        self.tau_correction_val = 0.009579
        #self.tau_correction_val = 0.0
        self.sorted_tauTheory = taus_of_funct(self.sorted_as, self.sorted_YTheory, self.tau_correction_val )
        self.sorted_tauTheory = [1.0 if tau > 1.0 else 0.0 if tau < 0.0 else tau for tau in self.sorted_tauTheory]
        self.sorted_tauMeas = taus_of_funct(self.sorted_as, self.sorted_YMeas, self.tau_correction_val )
        self.sorted_tauMeas = [1.0 if tau > 1.0 else 0.0 if tau < 0.0 else tau for tau in self.sorted_tauMeas]

        tau_scaling = 1.041
        self.sorted_tauMeas = [tau * tau_scaling for tau in self.sorted_tauMeas]

        Es_of_funct = lambda zs, OmM, OmL, Omk: np.sqrt(OmM * (1.0 + np.array(zs)) ** 3.0 + OmL + Omk * (1.0 + np.array(zs)) ** 2.0)
        self.sorted_Es = Es_of_funct(self.sorted_zs, OmM, OmL, Omk)

        self.canonicalTauOfzFunct = lambda zs,: [1.0 - integrate.quad(lambda zint: 1.0 / ((1.0 + zint) * Es_of_funct(zint, self.OmM, self.OmL, self.Omk)), 0 , z)[0] for z in zs]
        self.aSampling = np.linspace(0.0000, 1.0, n_t_interps)
        self.canon_tauSampling = [tau * tau_scaling for tau in self.canonicalTauOfzFunct([1.0 / a - 1.0 if a > 0.0 else np.inf for a in self.aSampling ]) ]

        self.tauToaInterp = interpolate.interp1d(self.canon_tauSampling, self.aSampling, bounds_error = False, fill_value = 0.0 )

        print('[min(self.sorted_tauMeas), max(self.sorted_tauMeas)] = ' + str([min(self.sorted_tauMeas), max(self.sorted_tauMeas)]))
        print('[min(self.canon_tauSampling), max(self.canon_tauSampling)] = ' + str([min(self.canon_tauSampling), max(self.canon_tauSampling)]))
        self.sorted_asMeasResiduals = np.array(self.sorted_as[1:] - self.tauToaInterp(self.sorted_tauMeas))

        self.n_bins = n_bins
        min_tau = min(self.sorted_tauMeas)
        max_tau = max(self.sorted_tauMeas)
        self.tau_bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
        self.tauMeas_bins, self.binned_tauMeas, self.binned_aMeas = bd.binData(self.sorted_tauMeas, self.sorted_as[1:], y_errs = None, bin_centers = None, bin_borders = self.tau_bin_edges, computed_value = 'mean', return_binned_xs = 1, trim = 0, print_frac = 1.1)
        self.binned_aMeas = self.binned_aMeas[0]
        self.tauMeas_bins, self.binned_tauMeas, self.binned_aMeasResiduals = bd.binData(self.sorted_tauMeas, self.sorted_asMeasResiduals, y_errs = None, bin_centers = None, bin_borders = self.tau_bin_edges, computed_value = 'mean', return_binned_xs = 1, trim = 0, print_frac = 1.1)
        self.binned_aMeasResiduals = self.binned_aMeasResiduals[0]
        #self.tauTheory_bins, self.binned_tauTheory, self.binned_aTheory = bd.binData(self.sorted_tauTheory, self.sorted_as[1:], y_errs = None, bin_centers = None, bin_borders = self.tau_bin_edges, computed_value = 'mean', return_binned_xs = 1, trim = 0, 0.05)
        #self.binned_aTheory = self.binned_aTheory[0]
        #self.mean_tauMeas = [np.mean(binned_taus) for binned_taus in self.binned_tauMeas]
        #self.mean_aMeas = [np.mean(binned_as) for binned_as in self.binned_aMeas]
        #self.mean_aMeasResiduals = [np.mean(binned_as) for binned_as in self.binned_aMeasResiduals]
        #self.mean_tauTheory = [np.mean(binned_taus) for binned_taus in self.binned_tauTheory]
        #self.mean_aTheory = [np.mean(binned_as) for binned_as in self.binned_aTheory]
        #self.ref_canonA = self.tauToaInterp(self.tauMeas_bins)
        #plt.scatter(self.tauMeas_bins, self.mean_aMeas)
        #plt.plot(self.tauMeas_bins, self.ref_canonA, c = 'k')
        #plt.show()

        #Now we take the derivatives of and smooth the data
        n_bins_for_deriv = 8
        delta_t = self.tauMeas_bins[1] - self.tauMeas_bins[0]
        self.binned_aMeasDerivs = [((self.binned_aMeasResiduals[min(n_bins-1, i + n_bins_for_deriv // 2 + 1)] - self.binned_aMeasResiduals[max(0, i - n_bins_for_deriv // 2)]) / ((min(n_bins-1, i + n_bins_for_deriv // 2 + 1) - max(0, i - n_bins_for_deriv // 2)) * delta_t) ) for i in range(n_bins)]

        #My attempts to exactly replicate the ksmooth command from MathCad from this link:
        # http://support.ptc.com/help/mathcad/en/index.html#page/PTC_Mathcad_Help/gaussian_kernel_smoothing.html
        gauss_kernel_funct = lambda ts: 1.0/ (np.sqrt(2.0 * np.pi) * 0.37) * np.exp(- ts ** 2.0 / (2.0 * 0.37 ** 2.0))
        ksmooth_funct = lambda vx, vy, b: [np.nansum([vy[j] * gauss_kernel_funct((vx_elem - vx[j]) / b) for j in range(len(vy))]) / np.nansum([gauss_kernel_funct((vx_elem - vx[j]) / b) for j in range(len(vx))]) for vx_elem in vx]

        self.time_cut_indeces = [ np.argmin(np.abs(np.array(self.tauMeas_bins) - self.time_cuts[0])), np.argmin(np.abs(np.array(self.tauMeas_bins) - self.time_cuts[1])) ]

        self.high_pass_smooth_aDerivs = ksmooth_funct(self.tauMeas_bins, self.binned_aMeasDerivs, 0.05)
        self.low_pass_smooth_aDerivs = ksmooth_funct(self.tauMeas_bins, self.binned_aMeasDerivs, 0.13)
        self.smooth_aDerivs_diff = np.array(self.high_pass_smooth_aDerivs) - np.array(self.low_pass_smooth_aDerivs)
        self.quantifyOscillationStrength()
