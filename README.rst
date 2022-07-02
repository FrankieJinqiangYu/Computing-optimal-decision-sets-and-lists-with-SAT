Learning optimal decision lists with SAT
-------------------------------------------------------------

*dlsat* is the MaxSAT-based approach in our paper 
`Learning Optimal Decision Sets and Lists with SAT
<https://www.jair.org/index.php/jair/article/download/12719/26747/>`_
to minimizing the total number of literals used in the decision list. 

**Getting Started**
`pandas
<https://pandas.pydata.org/>`_

`PySAT
<https://github.com/pysathq/pysat/>`_



An example of using the script to compute a complete perfect decision list is as follows:



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
