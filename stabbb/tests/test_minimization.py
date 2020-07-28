from .. import stab_BB
import numpy as np
from scipy.optimize import minimize, rosen, rosen_der

tol = 1e-6


def test_gradDesc_adaptative():
    """ Test whether the bbstab implementation gives the same result as 
    scipy.optimize for the Rosenbrock function"""
    x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    bestXstab, _ = stab_BB(x0, rosen, rosen_der)
    res = minimize(rosen, x0, method="Nelder-Mead", tol=1e-6)
    assert np.all(
        abs(res.x - bestXstab < tol)
    ), "gradient descent is not optimizing for adaptative deltaUpdateStrategy"


def test_gradDesc_constant():
    """ Test whether the bbstab implementation gives the same result as 
    scipy.optimize for the Rosenbrock function"""
    x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    bestXstab, _ = stab_BB(
        x0, rosen, rosen_der, bb=1, deltaUpdateStrategy="constant", deltaInput=2
    )
    res = minimize(rosen, x0, method="Nelder-Mead", tol=1e-6)
    assert np.all(
        abs(res.x - bestXstab < tol)
    ), "gradient descent is not optimizing for constant deltaUpdateStrategy"