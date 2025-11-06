import copy
import cvxpy as cp

def apply_canonical_condition(_X):
    """_summary_

    Args:
        _X (n,n): Base-pair probabilities

    Returns:
        optimal_X: Binary matrix satisfied secondary structure conditions
    """
    n = _X.shape[0]
    solver_name = "SCIP"
    A = _X
    B = copy.deepcopy(A)
    B[B>0] = 1
    B[B<=0] = 0
    B=B.astype(bool)
    X = cp.Variable((n, n), boolean=True)
    X.value = B

    constraints1 = [cp.sum(X[:, i]) <= 1 for i in range(n)]
    constraints2 = [X == X.T, X >= 0]

    obj = cp.Maximize(cp.sum(cp.multiply(A, X)))

    constraints = constraints1 + constraints2
    prob = cp.Problem(obj, constraints)
    result = prob.solve(solver_name)
    optimal_X = X.value

    return optimal_X
