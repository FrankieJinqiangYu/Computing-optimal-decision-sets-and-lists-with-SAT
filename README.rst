Learning Optimal Decision Lists with SAT
=============================================================

*dlsat* is the MaxSAT-based approach introduced in the paper 
`Learning Optimal Decision Sets and Lists with SAT
<https://www.jair.org/index.php/jair/article/download/12719/26747/>`_
to minimizing the total number of literals used in the decision list. 

Getting Started
---------------

`pandas
<https://pandas.pydata.org/>`_

`PySAT
<https://github.com/pysathq/pysat/>`_


Usage
-----

Before using *dlsat*, make sure you have the following Python packages installed:

``optdl.py`` provides a number of command-line options. To see these options by running:

::

   $ optdl.py -h

An example of using the script to compute a complete perfect decision list is as follows:


::

   $ optdl.py -a perfect --mode complete -v <dataset.csv>


An example of learning a separated perfect decision list is as follows:

::

   $ dlsolve.py -a perfect --clsorder maj --mode sep -v <dataset.csv>
   
``--clsorder`` enables to tool to sort classes by the number of items in each order. By default, it sorts in the ascending order. 
``--clsdown`` is optional. Activating this option such that classes are sorted in the descending order.

*dlsat* is also able to compute sparse decision lists. The examples of learning decision lists are as follow:

To compute complete sparse decision lists:

::

   $ optdl.py --approx 1 -a sparse --lambda 0.005 --mode complete -v <dataset.csv> 

To compute separated sparse decision lists:

::

   $ optdl.py --approx 1 -a sparse --lambda 0.005  --clsorder maj --mode sep -v <dataset.csv> 
   
Here, the value of ``--lambda`` is the regularization cost parameter, which equals 0.005 *by default*. It indicates how much adding a literal/rule to the decision list costs with respect to the overall accuracy increase.


Tweaking the MaxSAT solver
--------------------------

As some of the developed algorithms apply MaxSAT solving, it is sometimes
important to get the best performance of the underlying MaxSAT solver. The
``optdl.py`` tool provides a number of command-line options to tweak the
internal heuristics of the award-winning MaxSAT solver RC2 used in *dlsat*:

-  ``-1`` - to detect AtMost1 constraints
-  ``-b`` - to apply Boolean lexicographic optimization (BLO)
-  ``-m`` - to apply heuristic minimization of unsatisfiable cores
-  ``-t`` - to *trim* unsatisfiable cores (at most) a given number of times
-  ``-x`` - to *exhaust* unsatisfiable cores

You may want to use any combination of these. Also, note that *none of them*
are enabled by default. The techniques enabled by these command-line
parameters are detailed in `the paper describing RC2
<https://alexeyignatiev.github.io/assets/pdf/imms-jsat19-preprint.pdf>`__.
Read it if you are interested.
