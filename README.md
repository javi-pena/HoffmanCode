This is a python implementation of an algorithm to compute Hoffman constants as described in the Math Programming article "New characterizations of Hoffman constants for systems of
linear constraints", by Pena, Vera, Zuluaga.

The code takes as input a matrix $A\in \mathbb{R}^{m\times n}$ and computes the Hoffman constant $H(A)$ defined as follows.  For $b\in \mathbb{R}^m$ let $P_A(b):={x\in \mathbb{R}^n: Ax\le b}$.  Define
$$
H(A) = \sup_{b \in A(\mathbb{R}^n) + \mathbb{R}^m_+\atop x \not \in P_A(b)}\frac{\|(Ax-b)_+\|}{\text{dist}(x,P_A(b))}
$$

The main files are in the python file "Hoffman.py"
The python notebook "Hoffman.ipynb" illustrates the code in some instances of Examples 1, 2, and 3 from the above article.
