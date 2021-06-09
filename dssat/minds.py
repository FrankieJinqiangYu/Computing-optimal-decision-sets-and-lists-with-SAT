#!/usr/bin/env python
#-*- coding:utf-8 -*-
##
## minds.py
##
##  Created on: Dec 3, 2017
##      Author: Alexey Ignatiev
##      E-mail: alexey.ignatiev@monash.edu
##

# print function as in Python3
#==============================================================================
from __future__ import print_function
from check import ConsistencyChecker
from cover import Coverer
from hcover import HybCoverer
from data import Data
from minds1 import MinDS1
from mp92 import MP92
from options import Options
import os
from primer import Primer
import resource
from sat import SAT
from satlits import SATLits
from satlits3 import SATLits3
from mxsat import MaxSAT
from mxsatlits3 import MaxSATLits3
from mxsatlits4 import MaxSATLits4
import six
import sys


#
#==============================================================================
def show_info():
    """
        Print info message.
    """

    print('c MinDS: miner of decision sets')
    print('c author(s): Alexey Ignatiev [email:alexey.ignatiev@monash.edu]')
    print('')


#
#==============================================================================
def do_prime_based(data, options):
    """
        Run the prime-based approach.
    """

    if not options.approach.startswith('ph'):
        # phase1: enumeration of prime implicants
        primer = Primer(data, options)
        primes = primer.compute()
        print('c1 # of primes: {0}'.format(sum([len(p) for p in six.itervalues(primes)])))
        print('c1 prime time: {0:.4f}'.format(primer.time))
        print('')

        # phase2: computing smallest cover
        coverer = Coverer(data, primes, options)
    else:
        if options.approach in ('phybrid', 'phg'):
            options.approach = 'pgreedy'
        else:
            options.approach = 'pbased'

        # phase1: enumeration of prime implicants
        primer = Primer(data, options)
        primes = primer.compute()
        print('c1 # of primes: {0}'.format(sum([len(p) for p in six.itervalues(primes)])))
        print('c1 prime time: {0:.4f}'.format(primer.time))
        print('')

        # hybrid phase2 (column generation)
        coverer = HybCoverer(data, primer, primes, options)

    covers = coverer.compute()

    if isinstance(coverer, HybCoverer):
        print('c2 total # of primes: {0}'.format(sum([len(p) for p in six.itervalues(primes)])))

    print('c2 cover size: {0}'.format(sum([len(p) for p in six.itervalues(covers)])))

    if options.weighted:
        print('c2 cover wght: {0}'.format(coverer.cost))

    print('c2 cover time: {0:.4f}'.format(coverer.time))

    # save result to a CSV file
    if options.rdump:
        data.dump_result(primes, covers)


#
#==============================================================================
def do_sat_based(data, options):
    """
        Run the SAT-based approach.
    """

    if options.approach == 'mp92':
        solver = MP92(data, options)
    elif options.approach == 'sat':
        solver = SAT(data, options)
    elif options.approach == 'satlits':
        solver = SATLits(data, options)
    elif options.approach == 'satlits3':
        solver = SATLits3(data, options)
    elif options.approach in ('mxsatlits3', 'mxsat3', 'm3'):
        solver = MaxSATLits3(data, options)
    elif options.approach in ('mxsatlits4', 'mxsat4', 'm4'):
        solver = MaxSATLits4(data, options)
    elif options.approach in ('mxsat', 'mxsat4', 'mx'):
        solver = MaxSAT(data, options)
    else:
        solver = MinDS1(data, options)

    covers = solver.compute()

    if options.verb:
        # checking if there are default rules
        # and selecting the best among them
        wghts = []
        for lb, premise in six.iteritems(covers):
            # if len(premise) == 1 and premise[0] == []:
            wghts.append(tuple([lb, sum(data.wghts[i] for i in solver.samps[lb])]))

        if wghts:
            lb = max(filter(lambda p: len(covers[p[0]]) != 1 or covers[p[0]] != [], wghts), key=lambda p: p[1])[0]
            print('c1 cover: true => {0}'.format(': '.join(data.fvmap.opp[lb])))
        elif options.default:
            lb = max(wghts, key=lambda p: p[1])[0]
            print('c1 cover: true => {0}'.format(': '.join(data.fvmap.opp[lb])))

    print('c2 cover size: {0}'.format(sum([len(p) for p in six.itervalues(covers)])))
    print('c2 cover wght: {0}'.format(solver.cost))

    if hasattr(solver, 'accy'):
        print('c2 accy filtr: {0:.2f}%'.format(solver.accy))
    if hasattr(solver, 'accy_tot'):
        print('c2 accy total: {0:.2f}%'.format(solver.accy_tot))

    print('c2 cover time: {0:.4f}'.format(solver.time))

#
#==============================================================================
if __name__ == '__main__':
    # parsing command-line options
    options = Options(sys.argv)

    # making output unbuffered
    # sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    # showing head
    show_info()

    # parsing data
    if options.files:
        data = Data(filename=options.files[0], mapfile=options.mapfile,
                separator=options.separator, ranges=options.ranges)
    else:
        data = Data(fpointer=sys.stdin, mapfile=options.mapfile,
                separator=options.separator)

    if options.verb:
        print('c0 # of samps: {0} ({1} weighted)'.format(sum(data.wghts), len(data.samps)))
        print('c0 # of feats: {0} ({1} binary)'.format(len(data.names) - 1, len(list(filter(lambda x: x > 0, data.fvmap.opp.keys()))) - len(data.feats[-1])))
        print('c0 # of labls: {0}'.format(len(data.feats[-1])))

        used_time = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        print('c0 parse time: {0:.4f}'.format(used_time))
        print('')

    if options.noccheck == False:
        # phase0: consistency check
        checker = ConsistencyChecker(data, options)
        if checker.status and checker.do() == False:
            checker.remove_inconsistent()
            if options.verb:
                print('c0 data set is inconsistent')
                print('c0 filtering out {0} samples ({1} left)'.format(data.samps_filt, len(data.samps)))
                print('c0 filtering out {0} weights ({1} left)'.format(data.wghts_filt, sum(data.wghts)))
                print('c0 check time: {0:.4f}'.format(checker.time))
                print('')

            if options.cdump:
                checker.dump_consistent()

        if checker.status == False:
            print('c0 not enough classes => classification makes no sense')
            sys.exit(1)

    if options.approach in ('pbased', 'pgreedy', 'phybrid', 'phg', 'phb'):
        do_prime_based(data, options)
    else:  # sat
        do_sat_based(data, options)

    if options.verb:
        total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_utime
        print('c3 total time: {0:.4f}'.format(total_time))
