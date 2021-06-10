Learning optimal decision sets and lists with SAT
=============================================================

Decision sets and decision lists are two examples of the most easily explainable machine learning models. Here, we define size as the total number of literals in these rulebased models as opposed to earlier work that concentrates on the number of rules. In this paper, we develop approaches to computing minimum-size *perfect* decision sets and decision lists, which are perfectly accurate on the training data, and minimal in size, making use of modern SAT solving technology. We also provide a new method for determining optimal sparse alternatives, which trade off size and accuracy.


Approach to learning optimal decision sets with SAT
-------------------------------------------------------------

*minds* is a Python toolkit, which can be used for computing minimum size decisions sets, i.e. unordered sets of *if-then* rules [1]_. The toolkit represents a set of pure Python modules, which can be used in a Python script in the standard way through the provided API. Additionally, it contains an executable ``mds.py``, which can be applied for constructing a smallest size decision set for a training dataset in `CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__.

.. [1] Here the size can be defined as the number of rules of the decision set, or the total number of literals used.

Our novel SAT- and MaxSAT-based approach to minimizing the total number of literals used in the target decision set is implemted in *minds*. An example of how this can be done using the ``mds.py`` script follows:

::

   $ mds.py -a satl -s glucose3 -v <dataset.csv>

Here, one can replace the argument value ``'satl'`` with values ``satls`` to split the computation process by classes, or with values ``mxsatl`` and ``mxsatls`` to achieve the result by exploiting MaxSAT solvers (instead of iterative SAT solving).

Sparse decision sets can be constructed by running:

::

   $ mds.py -a sparse --lambda <float> -s glucose3 -v <dataset.csv>

Here, the value of ``'--lambda'`` is the regularization cost parameter, which equals 0.005 *by default*. It indicates how much adding a literal/rule to the decision set costs with respect to the overall accuracy increase.

There are many other approaches to learning optimal decisoin sets, detailed in `minds <https://github.com/alexeyignatiev/minds>`__.


Approach to learning optimal decision lists with SAT
-------------------------------------------------------------

*dlsat* is our MaxSAT-based approach to minimizing the total number of literals used in the decision list. An example of using the script to compute a complete perfect decision list is as follows:

::

   $ dlsolve.py --maxsat <dataset.csv>
   
An example of learning a separated perfect decision list is as follows:

::

   $ dlsolve.py -o maj -a <str> --maxsat --sep <dataset.csv>
   
The value of ``'-a'`` is either *asc* or *desc*, where the class order in the decision list is ascending/ descending regarding the number of items in each class.

*dlsat* is also able to compute sparse decision lists. The examples of learning decision lists are as follow:

To compute complete sparse decision lists:

::

   $ dlsolve.py -l <float> --maxsat --sparse <dataset.csv>
   
To compute separated sparse decision lists:

::

   $ dlsolve.py -l <float> -o <str> -a <str> --maxsat --sparse --sep <dataset.csv>
   
Here, the value of ``'-l'`` is the regularization cost parameter called *lambda*, which equals 0.005 *by default*. It indicates how much adding a literal/rule to the decision list costs with respect to the overall accuracy increase. The value ``'-o'`` is one of *maj*, *accuy* and *cost*, which represents the class order in the target decision list is based on the number of items/ accuracy/ cost of each class. The value of ``'-a'`` remains the same as perfect models, which is either *asc* or *desc*, representing the class in the decision list is in the ascending/ descending order regarding the metrics.
