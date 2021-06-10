The approach to learning optimal decision lists with SAT
-------------------------------------------------------------

*dlsat* is our MaxSAT-based approach to minimizing the total number of literals used in the decision list. An example of using the script to compute a complete perfect decision list is as follows:

::

   $ dlsolve.py --maxsat <dataset.csv>
   
An example of learning a separated perfect decision list is as follows:

::

   $ dlsolve.py --maxsat --sep -o maj -a <str> <dataset.csv>
   
The value of ``'-a'`` is either *asc* or *desc*, where the class order in the decision list is ascending/ descending regarding the number of items in each class.

*dlsat* is also able to compute sparse decision lists. The examples of learning decision lists are as follow:

To compute complete sparse decision lists:

::

   $ dlsolve.py --maxsat -l <float> <dataset.csv>
   
To compute separated sparse decision lists:

::

   $ dlsolve.py --maxsat -l <float> -o <str> -a <str> <dataset.csv>
   
Here, the value of ``'-l'`` is the regularization cost parameter called *lambda*, which equals 0.005 *by default*. It indicates how much adding a literal/rule to the decision list costs with respect to the overall accuracy increase. The value ``'-o'`` is one of *maj*, *accuy* and *cost*, which represents the class order in the target decision list is based on the number of items/ accuracy/ cost of each class. The value of ``'-a'`` remains the same as perfect models, which is either *asc* or *desc*, representing the class in the decision list is in the ascending/ descending order regarding the metrics.
