
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

import matplotlib.gridspec as gridspec
import scipy.stats

# TODO: multiple plot types: line plot with dist, just dist, maybe just line?
# also have option for horizontal

class Plot_Relaxation:
    """
    Plotting class for per frame relaxation data.
    """
    # any class attributes?

    def __init__(self, data, plot_type, ):
        """
        Parameters
        ----------

        """ 

    def formatting(self):
        plt.rcParams['figure.figsize']= (10,6)
        plt.rcParams.update({'font.size': 18})
        plt.rcParams["font.family"]="Sans-serif"
        plt.rcParams['font.sans-serif'] = 'Verdana'
        plt.rcParams['mathtext.default'] = 'regular'
        plt.rcParams['axes.linewidth'] = 2.5
        plt.rcParams['xtick.major.size'] = 6
        plt.rcParams['xtick.major.width'] = 2.5
        plt.rcParams['xtick.minor.size'] = 2
        plt.rcParams['xtick.minor.width'] = 2
        plt.rcParams['ytick.major.size'] = 6
        plt.rcParams['ytick.major.width'] = 2.5

    def pre_processing(file, time_units=10**6):
        """
        Processes raw time series data to appropriate units of time.
        """
        data = np.genfromtxt(file)
        # time units should mostly be in ps: convert to us
        time = np.divide(data[:,0], time_units)
        return np.vstack([time, data[:,1]])

    def avg_and_stdev(data_list):
        """
        Returns the average and stdev of multiple timeseries datasets.
        """
        # only the y values, x axis is time.
        data = [i[1] for i in data_list]
        return np.average(data, axis=0), np.std(data, axis=0)


    def line_plot(time, data, ylim=(0,5), ax=None, stdev=None, 
                label=None, leg_cols=5, color=None, ylabel=None):
        """
        Parameters
        ----------
        time : array
            Timeseries values.
        data : array
            Dataset values.
        ylim : tuple
            2 item tuple to set custom y limits.
        stdev : array
            Used to generate errors for the line plot.
        label : str
            Label the line plot.
        leg_cols : int
            Number of columns in the legend.

        """
        if ax is None:
            fig, ax = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios' : [20, 5]})
        else:
            fig = plt.gca()

        # line plot
        ax[0].plot(time, data, linewidth=1, alpha=0.8, label=label, color=color)
        #ax[0].axvline(2, color="k", lw=2, ls="--")
        ax[0].set_xlabel("Time ($\mu$s)", labelpad=10, fontweight="bold")
        #ax[0].set_ylabel(r"RMSD ($\AA$)", labelpad=10, fontweight="bold")
        ax[0].set_ylabel(r"19F to C=O Distance ($\AA$)", labelpad=10, fontweight="bold")
        ax[0].set_ylabel(ylabel, labelpad=10, fontweight="bold")
        ax[0].set_ylim(ylim)
        #ax[0].set_xticks(np.arange(0, time[-1] + (time[-1] / 5), time[-1] / 5), minor=True)
        ax[0].grid(alpha=0.5)

        # secondary kde distribution plot
        grid = np.arange(ylim[0], ylim[1], .01, dtype=float)
        density = scipy.stats.gaussian_kde(data)(grid)
        ax[1].plot(density, grid, color=color)
        # TODO: maybe normalize density to 1 and then set xticks np.arange(0, 1, 0.2)
        #ax[1].set_xticks(np.arange(0, 1.2, 0.2))
        ax[1].set_xticks(np.arange(0, np.max(density) + 0.5, 0.5))
        ax[1].xaxis.set_ticklabels([])
        ax[1].set_xlabel("Distribution", labelpad=28, fontweight="bold")
        ax[1].grid(alpha=0.5)

        # optionally plot the stdev using fill_between
        if stdev is not None:
            ax[0].fill_between(time, np.add(data, stdev), np.subtract(data, stdev), alpha=0.3, color=color)
        if label:
            ax[0].legend(loc=8, frameon=False, ncol=leg_cols, bbox_to_anchor=(0.5, -0.38))

        #fig.tight_layout()
        #fig.savefig("figures/test.png", dpi=300, transparent=False)

    def plot_19F_r1_r2(R, type):
        """
        Parameters
        ----------
        R : str
            Can be 'R1' or 'R2'.
        type : str
            Can be 'overall' (avg and stdev) or 'singles' (all individual datasets).
        """
        cmap = cm.Dark2
        norm = Normalize(vmin=0, vmax=3)

        fig, ax = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios' : [20, 5]})
        exp_r1 = {"w4f":1.99, "w5f":1.19, "w6f":1.25, "w7f":1.20}
        exp_r2 = {"w4f":109.1, "w5f":64.8, "w6f":63.0, "w7f":109.6}

        if R == "R1":
            ylim = (0, 5)
            index = 1
            for num, r in enumerate(exp_r1.values()):
                ax[0].axhline(y=r, xmin=0, xmax=1, color=cmap(norm(num)), linestyle="--")
                ax[1].axhline(y=r, xmin=0, xmax=1, color=cmap(norm(num)), linestyle="--")
                
        elif R == "R2":
            ylim = (60, 140)
            index = 2
            for num, r in enumerate(exp_r2.values()):
                ax[0].axhline(y=r, xmin=0, xmax=1, color=cmap(norm(num)), linestyle="--")
                ax[1].axhline(y=r, xmin=0, xmax=1, color=cmap(norm(num)), linestyle="--")

        else:
            raise ValueError("'r' must be 'R1' or 'R2'.")

        for num, sys in enumerate(["w4f", "w5f", "w6f", "w7f"]):
            ### overall avg and stdev
            if type == "overall":
                # SYS plot
                data = [pre_processing(f"ipq/{sys}/v{i:02d}/1us_noion/19F_R1_R2_newddsum.dat", 
                        time_units=10**3, index=index) for i in range(0, 5)]
                avg, stdev = avg_and_stdev(data)
                line_plot(data[0][0], avg, stdev=stdev, ax=ax, ylim=ylim, label=sys.upper(), 
                        leg_cols=4, ylabel="Relaxation Rate $s^{-1}$ " + f"({R})", color=cmap(norm(num)))

            ### individual datasets
            elif type == "singles":
                fig, ax = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios' : [20, 5]})
                for i in range(0, 5):
                    # SYS plot
                    data = pre_processing(f"ipq/{sys}/v{i:02d}/1us_noion/19F_R1_R2.dat", 
                                        time_units=10**3, index=index)
                    line_plot(data[0], data[index], ax=ax, ylim=ylim, 
                            ylabel="Relaxation Rate $s^{-1}$" + f"({R})", 
                            color=cmap(norm(num))
                            )
                fig.suptitle(sys)
                #plt.show()

            else:
                raise ValueError("Type arg must be 'overall' or 'singles'")

        fig.tight_layout()
        #plt.show()
        fig.savefig(f"figures/19F_{R}_5us_all_avg_stdev_newddsum.png", dpi=300, transparent=True)
