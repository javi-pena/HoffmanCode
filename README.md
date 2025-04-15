This is a python implementation of an algorithm to compute Hoffman constants as described in Section 3.2 of the Math Programming article "New characterizations of Hoffman constants for systems of linear constraints", by Pena, Vera, Zuluaga.

The code takes as input a matrix $A\in \mathbb{R}^{m\times n}$ and computes the Hoffman constant $H(A)$ defined as follows.  For $b\in \mathbb{R}^m$ let $P_A(b):=\{x\in \mathbb{R}^n: Ax\le b\}$.  Define

$
H(A) = \sup_{b \in A(\mathbb{R}^n) + \mathbb{R}^m_+\atop x \not \in P_A(b)}\frac{\|(Ax-b)_+\|}{\text{dist}(x,P_A(b))}.
$

For ease of computation, we assume that $\mathbb{R}^n$ and $\mathbb{R}^m$ are endowed with the $\ell_1$-norm and $\ell_\infty$-norm respectively.


As described in Section 3.2 of the article, the code computes $H(A)$ as well as sets $\mathcal F, \mathcal I$ certifying the optimality of $H(A)$.  

The syntax for the code is 

`H, FF, II, numitns = hoffman(A,maxiter)`

The output `numitns` records the total number of main iterations.  This output is $-1$ when the maximum number of iterations `\maxiter` is reached.  In that case the value of `H` is only a lower bound on $H(A)$.

The main files are in the python file "Hoffman.py"
This notebook illustrates the code in some instances of Examples 1, 2, and 3 from the above article.
