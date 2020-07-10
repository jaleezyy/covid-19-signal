#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def write_coverage_plot(depth_file, coverage_plot_name):
        if not os.path.exists(depth_file):
            return

        coverage = []
        for line in open(depth_file):
            t = line.split('\t')
            assert len(t) == 3
            coverage.append(int(float(t[2].strip("\n"))))

        coverage = np.array(coverage)
        assert np.all(coverage >= 0)

        n = len(coverage)
        assert n >= 1

        chunk_size = 2500
        nchunks = (n + chunk_size -1) // chunk_size

        kwds = {'wspace':0, 'hspace':0.2, 'bottom':0.02, 'top':0.98 }
        fig, axarr = plt.subplots(nchunks, 1, sharex=True, gridspec_kw=kwds)

        fig.set_figwidth(8)
        fig.set_figheight(0.75 * nchunks)

        for i, ax in enumerate(axarr):
            lo = i*chunk_size
            hi = min(n, (i+1)*chunk_size)
            label = f'{lo}-{hi}'

            ax.fill_between(np.arange(hi-lo), coverage[lo:hi] + 0.1, 1)
            for level in [ 1.0e1, 1.0e2, 1.0e3]:
                ax.plot([0,hi-lo], [level,level], ls=':', color='black')

            ax.set_yscale('log')
            ax.set_ylim(1.0, 3.0e4)
            ax.text(0.01, 0.95, label, verticalalignment='top', transform=ax.transAxes, color='red')

        print(f"Writing {coverage_plot_name}")
        plt.savefig(coverage_plot_name)
        plt.clf()

if __name__ == '__main__':

        write_coverage_plot(sys.argv[1], sys.argv[2])

