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


from bisect import bisect
import random
import math
import scipy
from isocalc import periodic_table

print_log_messages = True
_max_list_len = 100000

def log(msg):
    if print_log_messages == True:
        print(msg)

def frange(min_val, max_val, step):
    """Returns an iterator of floats
    in the range [<min_val>,<max_val>] including
    <min_val> but excluding <max_val> with
    an interval of <step>.
    >>> from pisotope import iso_gen
    >>> lst = list(iso_gen.frange(10,20,1))
    >>> min(lst)
    10
    >>> max(lst)
    19
    >>> len(lst)
    10
    """
    val = min_val
    while val < max_val:
        yield val
        val += step

def gaussian_iter(xs, a, b, fwhm):
    """Returns an iterator for application
    of the gaussian function with parameters
    <a>, <b> and <fwhm> to the list <xs>.
    """
    c = fwhm / 2.35482
    a *= 1.0 / (c*math.sqrt(2*math.pi))
    for x in xs:
        yield a * math.e ** (-(x - b)**2/(2*c**2))

def gaussian(xs, a, b, fwhm):
    """Applies the gaussian with the given
    parameters to the passed in list.
    """
    return list(gaussian_iter(xs, a, b, fwhm))


class SumFormula(object):

    def __init__(self, from_string=None, add_macros=None):
        self.atom_count_tuples = []
        self.macros = [('ACN', 'C2H3N1'),#
            ('MeCN', 'C2H3N1'), ('THF', 'C4H8O1'), ('DMSO', 'C2H6S1O1')]
        if from_string:
            self.from_string(from_string)

    def subtract(self, other_sum_formula):
        for sub_symb, _ in other_sum_formula.atom_count_tuples:
            assert(sub_symb in [sym for sym,_ in self.atom_count_tuples])
        for sub_symb, sub_count in other_sum_formula.atom_count_tuples:
            for i,tup in enumerate(self.atom_count_tuples):
                symb, count = tup
                if symb == sub_symb:
                    assert(count - sub_count >= 0)
                    self.atom_count_tuples[i] = symb, count - sub_count
        return self

    def add(self, other_sum_formula):
        for add_symb, add_count in other_sum_formula.atom_count_tuples:
            found_sym = False
            for i,tup in enumerate(self.atom_count_tuples):
                symb, count = tup
                if symb == add_symb:
                    found_sym = True
                    self.atom_count_tuples[i] = symb, count + add_count
            if not found_sym:
                self.atom_count_tuples.append((add_symb, add_count))
        return self


    def from_string(self, from_string):
        from_string = self._process_macros(from_string)
        print('--> reduced macro', from_string)
        symb = ''
        num = ''
        mult_fac = 1
        mult_fac_until = len(from_string)
        skip_after_paren = 0
        for i in range(len(from_string)):
            if skip_after_paren > 0:
                skip_after_paren -= 1
                assert(skip_after_paren >= 0)
                continue
            if from_string[i] == '(':
                if symb != '' and num != '':
                    self.atom_count_tuples.append( (symb, int(num)*mult_fac) )
                symb = ''
                num = ''

                j = i+1
                while from_string[j] != ')':
                    j += 1

                mult_fac_until = j
                j += 1

                new_mult_fac = ''
                while j < len(from_string) and from_string[j].isdigit():
                    new_mult_fac = new_mult_fac + from_string[j]
                    j += 1
                mult_fac = int(new_mult_fac)
                continue

            if from_string[i] == ')':
                pass

            if i > mult_fac_until:
                mult_fac = 1
                mult_fac_until = len(from_string)

            if from_string[i].isdigit():
                num += from_string[i]
            else:
                if num != '':
                    self.atom_count_tuples.append( (symb, int(num)*mult_fac) )
                    symb = ''
                    num = ''
                if from_string[i] != ')':
                    symb += from_string[i]
                else:
                    skip_after_paren = 0
                    j = i+1
                    while j < len(from_string) and from_string[j].isdigit():
                        skip_after_paren += 1
                        j += 1
        try:
            self.atom_count_tuples.append( (symb, int(num)) )
        except:
            pass

        print(self.atom_count_tuples)


    def _process_macros(self, s):
        for m in self.macros:
            s = s.replace(m[0], m[1])
        return s

class IsotopeDistribution(object):

    def __init__(self):
        self.pse = periodic_table.PeriodicTable()

    def generate(self, sum_formula, prob_threshold=0.001, fwhm=0.1, pad_left=3,
                 pad_right=3, interpolate_grid=0.005):
        print('Simulating isotopic distribution ...')
        if isinstance(sum_formula, str):
            sum_formula = SumFormula(from_string=sum_formula)
        ds = self._generate_dir(sum_formula, prob_threshold=prob_threshold)

        print('Simulating gaussians...')


        xs = [d[3] for d in ds]
        ys = [d[2] for d in ds]


        x_min = min(xs) - pad_left
        x_max = max(xs) + pad_right
        plot_xs = list(frange(x_min, x_max, interpolate_grid))
        plot_ys = [0.0 for _ in plot_xs]
        for i,peak_x in enumerate(xs):
            print('peak_x',peak_x)
            b = peak_x
            a = ys[i]
            gauss_ys = gaussian(plot_xs, a, b, fwhm)
            for i,py in enumerate(gauss_ys):
                plot_ys[i] += py
        return ds, plot_xs, plot_ys





    def _generate_dir(self, sum_formula, prob_threshold=0.001):
        def binom(i,j):
            return scipy.special.binom(j,i)
        lst = []

        pse = self.pse
        for atom_type, num_atoms in sum_formula.atom_count_tuples:
            elt = pse.elements[atom_type]
            # Iterate over isotopes indexed by nominal masses
            for nom_mass, isotope in elt.isotopes.items():
                    prob = isotope.prob
                    abs_mass = isotope.mass
                    # each iso dist is made up of atom types, nominal masses,
                    # the probability and the mass of all atoms together.
                    lst.append(([atom_type], [nom_mass], prob, abs_mass))
                    items_to_append = []
                    for itm in lst:
                        for i in range(1, num_atoms+1):
                            items_to_append.append((itm[0]+[atom_type]*i, itm[1]+[nom_mass]*i, itm[2]*((prob)**i)*binom(num_atoms-i,num_atoms), itm[3]+abs_mass*i))
                    # prevent addition of very unlikely isotope distributions
                    items_to_append = list(filter(lambda itm: itm[2] > prob_threshold, items_to_append))
                    # prevent duplicates
                    lst = lst + list(filter(lambda itm: itm not in lst, items_to_append))

                    assert(len(lst) < _max_list_len)
            lst = list(filter(lambda itm: self._contains_num_atoms_of_type(itm[0], num_atoms, atom_type), lst))
        #print(lst)
        return list(self._filter_according_to_sum_formula(lst, sum_formula))

    def _filter_according_to_sum_formula(self, lst, sf):
        sf = sf.atom_count_tuples
        for itm in lst:
            if all([self._contains_num_atoms_of_type(itm[0], tup[1], tup[0]) for tup in sf]):
                yield itm


    def _contains_num_atoms_of_type(self, lst_of_atom_types, num_atoms, atom_type):
        predicate = len([True for itm in lst_of_atom_types if str(itm) == str(atom_type)]) == num_atoms
        if predicate:
            print('!!!')
        return predicate



class RandomIsotopeDistribution(object):

    def __init__(self):
        self.pse = periodic_table.PeriodicTable()

    def generate(self, sum_formula, max_iterations=10000, float_accuracy=5, fwhm=0.5, peaks=True, interpolate_grid=0.025, pad_left=2, pad_right=2):
        print('Simulating isotopic distribution...')
        mass_intensity_dir = self.generate_dir(sum_formula, max_iterations=max_iterations,
                                               float_accuracy=float_accuracy)
        print('Simulating gaussians...')
        if peaks:
            x_min = min(list(mass_intensity_dir.keys())) - pad_left
            x_max = max(list(mass_intensity_dir.keys())) + pad_right
            xs = list(frange(x_min, x_max, interpolate_grid))
            ys = [0.0 for _ in xs]
            for peak_x in mass_intensity_dir.keys():
                b = peak_x
                a = mass_intensity_dir[peak_x]
                p_ys = gaussian(xs, a, b, fwhm)
                for i,py in enumerate(p_ys):
                    ys[i] += py
            return xs, ys
        return mass_intensity_dir


    def generate_dir(self, sum_formula, max_iterations=10000, float_accuracy=3):
        mass_intensity_dir = {}
        for it in range(max_iterations):
            mass = round(self.random_walk(sum_formula), float_accuracy)
            if mass in mass_intensity_dir:
                mass_intensity_dir[mass] += 1
            else:
                mass_intensity_dir[mass] = 1
        return mass_intensity_dir

    def old_take_random_isotope(self, elt):
        isotopes = list(elt.isotopes.values())
        P = [isotope.prob for isotope in isotopes]

        cdf = [P[0]]
        for i in range(1, len(P)):
            cdf.append(cdf[-1] + P[i])

        random_ind = bisect(cdf,random())
        return isotopes[random_ind]

    def _shuffle_two_lists(self,a,b):
        import random

        c = list(zip(a, b))

        random.shuffle(c)

        a, b = zip(*c)
        return a,b


    def take_random_isotope(self, elt):
        def weighted_choice(choices):
            total = sum(w for c, w in choices)
            r = random.uniform(0, total)
            upto = 0
            for c, w in choices:
                if upto + w >= r:
                    return c
                upto += w
            assert False, "Shouldn't get here"
        isotopes = list(elt.isotopes.values())
        P = [isotope.prob for isotope in isotopes]
        isotopes,P = self._shuffle_two_lists(isotopes,P)
        return weighted_choice(list(zip(isotopes,P)))

    def random_walk(self, sum_formula):
        mass = 0.0
        elements = self.pse.elements
        take_random_isotope = self.take_random_isotope
        for atom, count in sum_formula.atom_count_tuples:
            for i in range(count):
                mass += take_random_isotope(elements[atom]).mass
        return mass
