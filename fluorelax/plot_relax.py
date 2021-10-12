
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

    def __init__(self, data, plot_type, data_from_file=False):
        """
        Parameters
        ----------
        data : ndarray
            Array of dataset to load. Cols : Frame, R1, R2.
        plot_type : str
            Plot output, options are 'plot', 'dist'.
        data_from_file : bool
            When True, load the data arg as a filepath.
        """ 
        if data_from_file is True:
            self.data = np.genfromtxt(data)
        else:
            self.data = data
        self.fig, self.ax = plt.subplots()
        self.plot_type = plot_type
        self.cmap = cm.Dark2 # TODO: add as arg?
        self.norm = Normalize(vmin=0, vmax=3)

    def pre_processing(self, time_units=10**6):
        """
        Processes raw time series data to appropriate units of time.

        Parameters
        ----------
        time_units : int
            Factor to divide by to get the appropriate timescale for the dataset.
        """
        # time units should mostly be in ps: convert to us
        time = np.divide(self.data[:,0], time_units)
        # stack the new time axis and all other columns of the dataset
        self.data = np.vstack([time, self.data[:,1:]])

    def formatting(self):
        # TODO: add if plot_type == "line" and "dist"
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

    def dist_plot(self):
        """
        Plot just the distribution of R1 or R2 values.
        """
        ax = self.ax
        # secondary kde distribution plot
        grid = np.arange(self.ylim[0], self.ylim[1], .01, dtype=float)
        density = scipy.stats.gaussian_kde(self.data[:, self.index])(grid)
        #ax.plot(grid, density, color=self.color)
        ax.plot(grid, density)
        # TODO: maybe normalize density to 1 and then set xticks np.arange(0, 1, 0.2)
        #ax.set_xticks(np.arange(0, 1.2, 0.2))
        #ax.set_xticks(np.arange(0, np.max(density) + 0.5, 0.5))
        ax.set_xlabel("Relaxation ($s^{-1}$)", labelpad=28, fontweight="bold")

        # Remove the non-bottom spines
        for kw in ("left", "right", "top"):
           ax.spines[kw].set_visible(False)
        ax.axes.yaxis.set_visible(False)
        # self.ax.xaxis.set_ticks_position("bottom")
        # self.ax.yaxis.set_ticks_position("none")

        self.fig.tight_layout()

    def hist_plot(self):
        """
        Plot a histogram of R1 or R2 values.
        """
        ax = self.ax
        # secondary kde distribution plot
        ax.hist(self.data[:, self.index])
        # TODO: maybe normalize density to 1 and then set xticks np.arange(0, 1, 0.2)
        #ax.set_xticks(np.arange(0, 1.2, 0.2))
        #ax.set_xticks(np.arange(0, np.max(density) + 0.5, 0.5))
        ax.set_xlabel("Relaxation Rate ($s^{-1}$)", labelpad=28, fontweight="bold")

        # Remove the non-bottom spines
        for kw in ("left", "right", "top"):
           ax.spines[kw].set_visible(False)
        ax.axes.yaxis.set_visible(False)
        # self.ax.xaxis.set_ticks_position("bottom")
        # self.ax.yaxis.set_ticks_position("none")

        self.fig.tight_layout()

    # def line_plot(self, time, data, ylim=(0,5), ax=None, stdev=None, 
    #             label=None, leg_cols=5, color=None, ylabel=None):
    #     """
    #     Parameters
    #     ----------
    #     time : array
    #         Timeseries values.
    #     data : array
    #         Dataset values.
    #     ylim : tuple
    #         2 item tuple to set custom y limits.
    #     stdev : array
    #         Used to generate errors for the line plot.
    #     label : str
    #         Label the line plot.
    #     leg_cols : int
    #         Number of columns in the legend.
    #     """
    #     if ax is None:
    #         fig, ax = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios' : [20, 5]})
    #     else:
    #         fig = plt.gca()

    #     # line plot
    #     ax[0].plot(time, data, linewidth=1, alpha=0.8, label=label, color=color)
    #     #ax[0].axvline(2, color="k", lw=2, ls="--")
    #     ax[0].set_xlabel("Time ($\mu$s)", labelpad=10, fontweight="bold")
    #     #ax[0].set_ylabel(r"RMSD ($\AA$)", labelpad=10, fontweight="bold")
    #     ax[0].set_ylabel(r"19F to C=O Distance ($\AA$)", labelpad=10, fontweight="bold")
    #     ax[0].set_ylabel(ylabel, labelpad=10, fontweight="bold")
    #     ax[0].set_ylim(ylim)
    #     #ax[0].set_xticks(np.arange(0, time[-1] + (time[-1] / 5), time[-1] / 5), minor=True)
    #     ax[0].grid(alpha=0.5)

    #     # secondary kde distribution plot
    #     grid = np.arange(ylim[0], ylim[1], .01, dtype=float)
    #     density = scipy.stats.gaussian_kde(data)(grid)
    #     ax[1].plot(density, grid, color=color)
    #     # TODO: maybe normalize density to 1 and then set xticks np.arange(0, 1, 0.2)
    #     #ax[1].set_xticks(np.arange(0, 1.2, 0.2))
    #     ax[1].set_xticks(np.arange(0, np.max(density) + 0.5, 0.5))
    #     ax[1].xaxis.set_ticklabels([])
    #     ax[1].set_xlabel("Distribution", labelpad=28, fontweight="bold")
    #     ax[1].grid(alpha=0.5)

    #     # optionally plot the stdev using fill_between
    #     if stdev is not None:
    #         ax[0].fill_between(time, np.add(data, stdev), np.subtract(data, stdev), alpha=0.3, color=color)
    #     if label:
    #         ax[0].legend(loc=8, frameon=False, ncol=leg_cols, bbox_to_anchor=(0.5, -0.38))

    #     #fig.tight_layout()
    #     #fig.savefig("figures/test.png", dpi=300, transparent=False)

    def plot_r1(self):
        """
        Main method for plotting the R1 values.
        """
        exp_r1 = {"w4f":1.99, "w5f":1.19, "w6f":1.25, "w7f":1.20}
        self.ylim = (0, 5)
        self.index = 1
        for num, r in enumerate(exp_r1.values()):
            self.ax.vlines(x=r, ymin=0, ymax=1, color=self.cmap(self.norm(num)), linestyle="--")
        self.dist_plot()

    def plot_r2(self):
        """
        Main method for plotting the R2 values.
        """
        self.ylim = (60, 140)
        self.index = 2
        exp_r2 = {"w4f":109.1, "w5f":64.8, "w6f":63.0, "w7f":109.6}
        for num, r in enumerate(exp_r2.values()):
            self.ax.vlines(x=r, ymin=0, ymax=1, color=self.cmap(self.norm(num)), linestyle="--")
        self.dist_plot()


    # split into 2 methods: r1 and r2 plots
    # def plot_19F_r1_r2(self, R, type):
    #     """
    #     Parameters
    #     ----------
    #     R : str
    #         Can be 'R1' or 'R2'.
    #     type : str
    #         Can be 'overall' (avg and stdev) or 'singles' (all individual datasets).
    #     """

    #     for num, sys in enumerate(["w4f", "w5f", "w6f", "w7f"]):
    #         ### overall avg and stdev
    #         if type == "overall":
    #             # SYS plot
    #             data = [pre_processing(f"ipq/{sys}/v{i:02d}/1us_noion/19F_R1_R2_newddsum.dat", 
    #                     time_units=10**3, index=index) for i in range(0, 5)]
    #             avg, stdev = avg_and_stdev(data)
    #             line_plot(data[0][0], avg, stdev=stdev, ax=ax, ylim=ylim, label=sys.upper(), 
    #                     leg_cols=4, ylabel="Relaxation Rate $s^{-1}$ " + f"({R})", color=cmap(norm(num)))

    #         ### individual datasets
    #         elif type == "singles":
    #             fig, ax = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios' : [20, 5]})
    #             for i in range(0, 5):
    #                 # SYS plot
    #                 data = pre_processing(f"ipq/{sys}/v{i:02d}/1us_noion/19F_R1_R2.dat", 
    #                                     time_units=10**3, index=index)
    #                 line_plot(data[0], data[index], ax=ax, ylim=ylim, 
    #                         ylabel="Relaxation Rate $s^{-1}$" + f"({R})", 
    #                         color=cmap(norm(num))
    #                         )
    #             fig.suptitle(sys)
    #             #plt.show()

    #         else:
    #             raise ValueError("Type arg must be 'overall' or 'singles'")

    #     fig.tight_layout()
    #     #plt.show()
    #     fig.savefig(f"figures/19F_{R}_5us_all_avg_stdev_newddsum.png", dpi=300, transparent=True)
