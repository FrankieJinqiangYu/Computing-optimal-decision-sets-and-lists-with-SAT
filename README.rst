Learning optimal decision sets and lists with SAT
=============================================================

Decision sets and decision lists are two examples of the most easily explainable machine learning models. Here, we define size as the total number of literals in these rulebased models as opposed to earlier work that concentrates on the number of rules. In this paper, we develop approaches to computing minimum-size *perfect* decision sets and decision lists, which are perfectly accurate on the training data, and minimal in size, making use of modern SAT solving technology. We also provide a new method for determining optimal sparse alternatives, which trade off size and accuracy.

*minds* is a Python toolkit, which can be used for computing minimum size decisions sets, i.e. unordered sets of *if-then* rules [1]_. The toolkit represents a set of pure Python modules, which can be used in a Python script in the standard way through the provided API. Additionally, it contains an executable ``mds.py``, which can be applied for constructing a smallest size decision set for a training dataset in `CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__.

.. [1] Here the size can be defined as the number of rules of the decision set, or the total number of literals used.

Our `recent CP'20 paper<https://alexeyignatiev.github.io/assets/pdf/yisb-cp20-preprint.pdf>`__ proposed a novel SAT- and MaxSAT-based approach to minimizing the total numberof literals used in the target decision set. An example of how this can be done using the ``mds.py`` script follows:

::

   $ mds.py -a satl -s glucose3 -v <dataset.csv>

Here, one can replace the argument value ``'satl'`` with values ``satls`` to split the computation process by classes, or with values ``mxsatl`` and ``mxsatls`` to achieve the result by exploiting MaxSAT solvers (instead of iterative SAT solving).

Sparse decision sets can be constructed by running:

::

   $ mds.py -a sparse --lambda <float> -s glucose3 -v <dataset.csv>

Here, the value of ``'--lambda'`` is the regularization cost parameter, whichequals 0.005 *by default*. It indicates how much adding a literal/rule to the
decision set costs with respect to the overall accuracy increase (see `the paper<https://alexeyignatiev.github.io/assets/pdf/yisb-cp20-preprint.pdf>`__ for details).
