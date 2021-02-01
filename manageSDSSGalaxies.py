import cantrips as c

class SDSSDataManager:

    def generateSDSSHistograms(self, fields = [0, 2, 3, 4, 5, 6, 7, 8, 9], hist_cropped_field = 1, hist_full_field = 1, plot_array_rows = 5,
                                 hist_colors = ['k','r'], ra_index = 1, dec_index = 2, redshift_index = 11, n_bins = 50, figsize = (10,8) ):
        sdss_loads = {field:c.readInColumnsToList(self.SDSSPhotozs[field],self.SDSSDir, n_ignore = 2, delimiter= ',') for field in fields}
        f, axarr = plt.subplots(plot_array_rows, (len(fields) // plot_array_rows) + 1, sharex = True, sharey = True, figsize = figsize)
        plt.subplots_adjust(hspace=0.0)
        plt.suptitle('SDSS Photo-zs in PS1MD Fields (white = 7 deg. diam. circle around field; orange = just field)')
        for i in range(len(fields)):
            ax = axarr[plot_array_rows - i % plot_array_rows - 1, i // plot_array_rows]
            field = fields[i]
            ras = [float(ra) for ra in sdss_loads[field][ra_index]]
            decs = [float(dec) for dec in sdss_loads[field][dec_index]]
            zs = [float(z) for z in sdss_loads[field][redshift_index]]
            if hist_full_field:
                ax.hist(zs, edgecolor = hist_colors[0], bins = n_bins, fill = False)
            if hist_cropped_field:
                field_bounds = self.fields[field]
                ax.hist([zs[i] for i in range(len(zs)) if (ras[i] > field_bounds[0] and ras[i] < field_bounds[1]
                                                           and decs[i] > field_bounds[2] and decs[i] < field_bounds[3]) ],
                        edgecolor = hist_colors[0], bins = n_bins, fill = hist_colors[1])
            if i % plot_array_rows == 0: ax.set_xlabel('SDSS photometric redshift')
            if i // plot_array_rows == 0: ax.set_ylabel('# Galaxies in Bin')
            ax.text(0.5, 700, 'PS1MD Field ' + str(field))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        return 1

    def __init__(self):
        print ('SDSS galaxy plotter initiated. ')
