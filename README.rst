# Computing-optimal-decision-sets-and-lists-with-SAT

*minds* is a Python toolkit, which can be used for computing minimum size
decisions sets, i.e. unordered sets of *if-then* rules [1]_. The toolkit
represents a set of pure Python modules, which can be used in a Python script
in the standard way through the provided API. Additionally, it contains an
executable ``mds.py``, which can be applied for constructing a smallest size
decision set for a training dataset in `CSV
<https://en.wikipedia.org/wiki/Comma-separated_values>`__.

.. [1] Here the size can be defined as the number of rules of the decision
   set, or the total number of literals used.
