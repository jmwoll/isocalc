# This file is licensed under the BSD 3-Clause License.
#
# Copyright (c) 2018, Jan Wollschläger
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from isocalc.iso_gen import RandomIsotopeDistribution, SumFormula
from matplotlib import pyplot as plt


def plot_isotope_distribution(s,save_as=None,charge=None,start=None,end=None,fragment_loss=None,adduct=None,retrieve_data=False,max_iterations=10000, float_accuracy=5, fwhm=0.1, peaks=True, interpolate_grid=0.025, pad_left=2, pad_right=2,show=True):
    """ Plots the isotope distribution for the string <s>
    with the charge state <charge> in the range <start>,<end>.
    Plot can be saved by specifying <save_as>
    """
    xy = isotope_distribution(s,fragment_loss=fragment_loss,adduct=adduct,max_iterations=max_iterations, float_accuracy=float_accuracy, fwhm=fwhm, peaks=peaks, interpolate_grid=interpolate_grid, pad_left=pad_left, pad_right=pad_right)
    if charge is not None and int(abs(charge)) != 1:
        abs_charge = abs(charge)
        print('z = ',charge)
        x_lst = xy[0]
        xy = ([itm/abs_charge for itm in x_lst],xy[1])
    plot_mass_spectrum(xy,direct_data_feed=True,save_as=save_as,show=show)
    if retrieve_data:
        return xy


_isotope_distribution_memo = {}
def isotope_distribution(s,fragment_loss=None,adduct=None,max_iterations=10000, float_accuracy=5, fwhm=0.1, peaks=True, interpolate_grid=0.025, pad_left=2, pad_right=2, memoize=None):
    """Generates the isotope distribution for
    the string <s> and returns <ds>, <plot_xs>, <plot_ys>
    where <ds> is a dictionary containing the key value
    pairs while plot_xs and plot_ys are the xy data of
    a simulation of the isotope distribution.
    """
    global _isotope_distribution_memo
    if memoize is None:
        memoize = True
    if memoize:
        if s in _isotope_distribution_memo:
            return _isotope_distribution_memo[s]
    if isinstance(s,str):
        s = SumFormula(s)
    if fragment_loss is not None:
        if isinstance(fragment_loss,str):
            fragment_loss = SumFormula(fragment_loss)
        s = s.subtract(fragment_loss)
    if adduct is not None:
        if isinstance(adduct, str):
            adduct = SumFormula(adduct)
        s = s.add(adduct)

    iso_dist = RandomIsotopeDistribution()
    rslt = iso_dist.generate(s,max_iterations=max_iterations, float_accuracy=float_accuracy, fwhm=fwhm, peaks=peaks, interpolate_grid=interpolate_grid, pad_left=pad_left, pad_right=pad_right)
    if memoize:
        _isotope_distribution_memo[s] = rslt
    return rslt




def scale(ys, scale=1):
    return list(map(lambda y: y*scale, ys))


def norm(xs):
    lstmax, lstmin = max(xs), min(xs)
    return list(map(lambda x: (x - lstmin) / float(lstmax - lstmin) ,xs))


def annotate_ms(xs, ys, int_thresh=.5, margin_x=0.05, margin_y=0.05, decimal_places=1):
    # annotate the given mass spectrum, such that the peaks
    # are annotated by their respective m/z
    max_y = max(ys)
    format_str = '{0:.'+str(decimal_places)+'f}' if decimal_places is not None else '{}'
    def annotate_peak(x,y):
        #plt.gca().annotate('{}'.format(x),(x,y),textcoords='data',horizontal_alignment='center')
        plt.gca().text(x,y+margin_y*max_y,format_str.format(x),horizontalalignment="center")

    thresh = int_thresh * max_y
    visited = {}
    for i,x in enumerate(xs):
        y = ys[i]
        try:
            if y >= max(ys[i-10:i]+ys[i+1:i+10]):
                if not(x in visited) and y > thresh:
                    visited[x] = True
                    #print('annotate {} {}'.format(x,y))
                    annotate_peak(x,y)
        except IndexError:
            pass

_plot_mass_spectrum_cache, _plot_mass_spectrum_cache_size = {}, 4
def plot_mass_spectrum(sample, start=None, end=None, int_thresh=.5, margin_x=0.05, margin_y=0.05,
                       save_as=None, scale_relative=True, decimal_places=1, remove_frequencies=False,
                       title=None, add_background_frequencies=False,
                       report=False, process_ys=None, process_xs=None,
                       process_xy=None, direct_data_feed=False, report_n=None,
                       get_data=None, integrate=None, show=None, rasterize=None,
                       fig_size=None, cache=None, annotate=None):
    global _plot_mass_spectrum_cache
    if len(_plot_mass_spectrum_cache) > _plot_mass_spectrum_cache_size:
        _plot_mass_spectrum_cache = {}

    xs,ys = None,None
    if direct_data_feed:
        assert(sample is not None)
        xs,ys = sample
    else:
        if cache is None:
            cache = True
        if cache and sample in _plot_mass_spectrum_cache:
            xs, ys = _plot_mass_spectrum_cache[sample]
        else:
            try:
                xs, ys = load_ms(sample, sep=',')
            except:
                xs, ys = load_ms(sample, sep='\t')
            if cache:
                _plot_mass_spectrum_cache[sample] = (xs,ys)

    assert(xs != None and ys != None)
    org_xs, org_ys = [itm for itm in xs], [itm for itm in ys]
    if process_xy is not None:
        assert(process_xs is None and process_ys is None)
        xs, ys = process_xy(xs, ys)
    elif process_xs is not None:
        xs = list(process_xs(xs))
    elif process_ys is not None:
        ys = list(process_ys(ys))

    if remove_frequencies:
        xs, ys = remove_background_frequencies(xs, ys,
                                add_background_freqs=add_background_frequencies)

    if report:
        show_report(xs, ys, report_n=report_n)

    if fig_size is not None:
        #plt.figure(figsize=(8*1.5,6*1.5))
        plt.figure(figsize=fig_size)

    if (start is not None) and (end is not None):
        xs, ys = cut_xy(xs, ys, start, end)

    if scale_relative:
        ys = norm(ys)
        ys = scale(ys,scale=100)

        plt.ylim([0,105])

    if start is None:
        start = min(xs)
    if end is None:
        end = max(xs)
    ax = plt.gca()
    ax.get_yaxis().set_tick_params(right=False, which='both', direction='out')
    ax.get_xaxis().set_tick_params(top=False, which='both', direction='out')
    plt.xlabel('m / z', style='italic')
    if scale_relative:
        plt.ylabel('rel. intensity / %')
    else:
        plt.ylabel('arb. intensity / counts')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if title is not None:
        plt.title(title)
    #xs,ys = cut_xy(xs,ys,start,end)
    if rasterize is not None:
        if rasterize:
            rasterize = True
    else:
        rasterize = False
    if rasterize:
        plt.plot(xs,ys,color='black',rasterized=rasterize)
    else:
        plt.plot(xs,ys,color='black')
    if annotate is None:
        annotate = True
    if show is None:
        show = True
    if show and annotate:
        annotate_ms(xs,ys,int_thresh=int_thresh,margin_x=margin_x,margin_y=margin_y,decimal_places=decimal_places)

    if save_as is not None:
        plt.savefig(save_as)

    if show == False:
        plt.clf()
    else:
        plt.show()

    data = {}
    if integrate is not None:
        assert(len(integrate) == 2)
        int_from = integrate[0]
        int_to = integrate[1]
        assert(int_from is not None)
        assert(int_to is not None)
        intxs, intys = cut_xy(org_xs, org_ys, int_from, int_to)
        data['int'] = np.trapz(intys, x=intxs)
        get_data = True

    if get_data is not None:
        if get_data:
            data['xs'] = xs
            data['ys'] = ys
            return data
