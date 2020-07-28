import numpy as np
import pdb

def stab_BB(
    x0,
    costFn,
    gradFn,
    bb = 1,
    deltaUpdateStrategy = 'adaptative',
    deltaInput = 1e6,
    c = 0.2,
    maxIt = 10000,
    tol = 1e-7,
    verbose=False,
):
    """ Stabilized Barzilai-Borwein Method for unconstrained optimization

    Implementes the algorithm described in [1] the paper Burdakov, Oleg & Dai, 
    Yu-Hong & Huang, Na. (2019). Stabilized Barzilai-Borwein.

    Parameters
    ----------
    x0 : np.array, shape(n,)
        First initial point
    costFn: callable
        The objective function to be minimized.

            ``fun(x) -> float``

        where ``x`` is an 1-D array with shape (n,)
    gradFn : callabe
        Function to evaluate the gradient of the cost function being minimized
    bb: {1,2}, optional (default 'adaptative')
        Whether to use BB1 or BB2 method
    deltaUpdateStrategy: {'adaptative','constant'}
        Whether to use a constant or adaptative delta parameter
    deltaInput: number, optional (default 1e6)
        If ``deltaUpdateStrategy=='adaptative'`` it's the delta used until the
        fourth iteration. For this case it's suggested to use a large value like
        1e6.
        Else it will be the delta used throughout the experiment. 
        The paper was tested setting this value to be 2 or 0.01, 0.1 or 1.0.
    c: number, optional (default 1e6)
        Paramater for adaptative choice of delta. Must be set if
        ``deltaUpdateStrategy=='adaptative'``.
        The paper was tested using 0.1,0.2,0.25 or 0.3 for quadratic cost 
        functions and also using 0.1, 0.5 or 1.0 for nonquadratic ones.
    maxIt : int, optional (defalt 500)
        Maximum number of iteration
    tol : float, optional (default 1e-7)
        Value that will cause the algorithm to stop
    verbose : bool, optional (defalt False)
        Whether to print partial results after each iteration

    Returns
    -------
    bestX, (xHistory, alphaHistory, normGradHistory)
        ``bestX`` is (n,) ndarray representing the best evaluated point
        ``xHistory`` is a list of all evaluated points
        ``alphaHistory`` is float list of all evaluate alphas
        ``normGradHistory`` is float list of all evaluated norm of gradients

    References
    ----------
    .. [1] Burdakov, Oleg & Dai, Yu-Hong & Huang, Na. (2019). Stabilized 
        Barzilai-Borwein.

    Examples
    ---------
    Let's use the Rosenbrock function.

    >> from scipy.optimize import rosen, rosen_der
    >> from stabbb import stab_BB
    >> x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    >> bestX, _ = stab_BB(x0, rosen, rosen_der)
    >> bestX
    array([ 1.,  1.,  1.,  1.,  1.])
    """

    class StabBBInt:
        """ Internals of the function 
        The motivation for using this internal class is to be able to return a
        value even if the user makes a keyboard interruption.
        """

        def __init__(self):
            #: History of evaluated x's
            self.x = [x0, self.backtracking(x0, gradFn(x0))]
            #: History of evaluated alphas
            self.alpha = []

        def backtracking(self, x0, g0):
            """ Suggested algorithm to avoid poor choices of 2 initial points """
            alpha0 = 1 / np.linalg.norm(x0, np.inf)
            s0 = -alpha0 * g0
            x1 = x0 + s0
            while costFn(x1) > costFn(x0):
                s0 = s0 / 4
                x1 = x0 + s0
            return x1

        def mainAlg(self):
            delta = deltaInput
            #: g_{k-1} for the algorithm
            gkant = gradFn(x0)
            #: g_{k} for the algorithm
            gk = gradFn(self.x[1])
            #: History of gradient norm
            self.normGrad = [np.linalg.norm(gk)]
            for k in range(int(maxIt)):
                if self.normGrad[-1] > tol:
                    #: $s_{k-1} \leftarrow x_{k-1} - x_{k-2}$
                    skant = self.x[-1] - self.x[-2]
                    #: $y_{k-1} \leftarrow  g_k - g_{k-1}
                    ykant = gk - gkant
                    
                    # Compute $\alpha_k$
                    #
                    alpha_bb = self.calcAlpha(skant, ykant)
                    # Adaptation of method to deal with nonconvex functions
                    alpha_bb_corrected = self.correctAlpha(alpha_bb, skant, ykant)
                    # Applies the proposed stabilization
                    alphak = min([alpha_bb_corrected, delta / np.linalg.norm(gk),])

                    #: Set $x_{k+1} \leftarrow x_k - \alpha_k g_k$
                    xk = self.x[-1] - alphak * gk
                    
                    # Update $g$ 
                    gkant = gk
                    gk = gradFn(xk)

                    # Update $\Delta$ depending on strategy
                    if deltaUpdateStrategy == "adaptative":
                        delta = self.deltaUpdateFn(delta, k)

                    # Store processed states
                    self.alpha.append(alphak)
                    self.x.append(xk)
                    self.normGrad.append(np.linalg.norm(gk))

                    if verbose:
                        print(f"New norm of gradient is {self.normGrad[-1]}")

        def calcAlpha(self, skant, ykant):
            if bb == 1:
                return (skant.transpose() @ ykant) / (ykant.transpose() @ ykant)
            elif bb == 2:
                return (skant.transpose() @ skant) / (skant.transpose() @ ykant)
            else:
                raise Exception("parameter bb is not within range")

        def correctAlpha(self, alpha_bb, skant, ykant):
            if alpha_bb <= 0:
                # eq 4.2
                # Although the paper says this correction for nonconvex
                # problems depends on sk for the current iteration,
                # it's not available up to this point. Will use k-1
                # regardless
                return np.linalg.norm(skant) / np.linalg.norm(ykant)
            else:
                return alpha_bb

        def deltaUpdateFn(self, lastDelta, k):
            if k == 4:
                # Adaptative choice of parameter delta as equation (4.1)
                return c * min(
                    [
                        np.linalg.norm(self.x[4] - self.x[3]),
                        np.linalg.norm(self.x[3] - self.x[2]),
                        np.linalg.norm(self.x[2] - self.x[1]),
                    ]
                )
            else:
                return lastDelta

    stabint = StabBBInt()

    try:
        stabint.mainAlg()
    except KeyboardInterrupt:
        pass
    bestX = stabint.x[np.argmin(stabint.normGrad) + 1]
    xHistory, alphaHistory, normGradHistory = (
        stabint.x,
        stabint.alpha,
        stabint.normGrad,
    )
    return bestX, (xHistory, alphaHistory, normGradHistory)
