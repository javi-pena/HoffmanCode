"""
Created on Sat Feb  4 09:41:14 2023

@author: Javier Peña
"""

import cvxpy as cp
import numpy as np

def test(A):
    """
    Solves the convex optimization problem:
    
        minimize t
        subject to:
        A * x <= t
        A * x >= -t
        sum(x) = 1
        x >= 0, t >= 0

    Parameters:
        A (numpy.ndarray): Input matrix.

    Returns:
        x (numpy.ndarray): Optimal solution vector.
        t (float): Optimal value of t.
    """
    n = A.shape[1]  # Number of variables
    m = A.shape[0]  # Number of constraints

    # Define variables
    x = cp.Variable(n)
    t = cp.Variable()

    # Define constraints
    constraints = [
        A @ x <= t * np.ones(m),
        A @ x >= -t * np.ones(m),
        cp.sum(x) == 1,
        t >= 0,
        x >= 0
    ]

    # Define optimization problem
    problem = cp.Problem(cp.Minimize(t), constraints)

    # Solve using Gurobi
    problem.solve(solver=cp.GUROBI)

    return x.value, t.value

def findJ(FF, II, oldJ):
    m = max(FF.shape[1] if FF.size else 0, II.shape[1] if II.size else 0)
    
    x = cp.Variable(m, boolean=True)  # Binary variables

    # Objective: Maximize sum(x)
    objective = cp.Maximize(cp.sum(x))

    # Constraints
    constraints = []
        
    if FF.shape[0] == 0:
        constraints.append(II @ x <= II @ np.ones(m) - 1)
    elif II.shape[0] == 0:
        constraints.append((FF - 1) @ x <= -np.ones(FF.shape[0]))
    else:
        constraints.append((FF - 1) @ x <= -np.ones(FF.shape[0]))
        constraints.append(II @ x <= II @ np.ones(m) - 1)

    # Extra constraint
    constraints.append(cp.sum(x) <= np.sum(oldJ))

    # Solve using Gurobi
    problem = cp.Problem(objective, constraints)
    problem.solve(solver=cp.GUROBI)
    try:
        return x.value.round().astype(int)  # Convert to integer array
    except:
        return np.array([])

def hoffman(A, maxiter):
    """
    Compute Hoffman constant H for matrix A.
    Also computes FF, II providing certificates of surjectivity and non-surjectivity.
    
    Parameters:
        A (numpy.ndarray): Input matrix.
        maxiter (int): Maximum number of iterations.
    
    Returns:
        H (float): Hoffman constant.
        FF (numpy.ndarray): Collection of surjectivity certificates.
        II (numpy.ndarray): Collection of non-surjectivity certificates.
        iternumber (int): Number of iterations performed.
    """
    m = A.shape[0]
    FF = np.empty((0, m), dtype=int)
    II = np.empty((0, m), dtype=int)
    J = np.ones(m, dtype=int)
    H = 0
    iternumber = 0

    while J.size > 0 and iternumber < maxiter:
        iternumber += 1
        AA = A[J > 0, :].T  # Select rows where J > 0 and transpose
        
        y, t = test(AA)

        if t > 0:
            # J is A-surjective, add to FF and update H
            FF = np.vstack([FF, J])
            H = max(H, 1 / t)
            
            if 1 / t >= H:
                bestJ = J
                besty = y
        else:
            # Get a certificate of non-surjectivity for J and add it to II
            yy = np.zeros(m)
            yy[J > 0] = y
            
            # I is a binary vector indicating positions where yy > 0
            I = np.zeros(m, dtype=int)
            I[yy > 0] = 1
            II = np.vstack([II, I])
        
        J = findJ(FF, II, J)
    
    if iternumber == maxiter:
        iternumber = -1
    
    return H, FF, II, iternumber
