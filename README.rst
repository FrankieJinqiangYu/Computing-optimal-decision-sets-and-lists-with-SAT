Learning optimal decision sets and lists with SAT
=============================================================

Decision sets and decision lists are two examples of the most easily explainable machine learning models. Here, we define size as the total number of literals in these rulebased models as opposed to earlier work that concentrates on the number of rules. In this paper, we develop approaches to computing minimum-size *perfect* decision sets and decision lists, which are perfectly accurate on the training data, and minimal in size, making use of modern SAT solving technology. We also provide a new method for determining optimal sparse alternatives, which trade off size and accuracy.

*minds* is a Python toolkit, which can be used for computing minimum size decisions sets, i.e. unordered sets of *if-then* rules [1]_. The toolkit represents a set of pure Python modules, which can be used in a Python script in the standard way through the provided API. 

.. [1] Here the size can be defined as the number of rules of the decision
   set, or the total number of literals used.
