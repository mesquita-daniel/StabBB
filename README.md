# Stabilized Barzilai-Borwein method

Implements gradient descent using Stabilized Barzilai-Borwein method.

It's based on the following paper: [Burdakov, Oleg & Dai, Yu-Hong & Huang, Na. (2019). Stabilized Barzilai-Borwein](https://arxiv.org/pdf/1907.06409.pdf).

The Barzilai-Borwein (BB) method is a popular and efficient tool for solving large-scale unconstrained optimization problems. Its search direction is the same as for the steepest descent (Cauchy) method, but its stepsize rule is different. Owing to this, it converges much faster than the Cauchy method. A feature of the BB method is that it may generate too long steps, which throw the iterates too far away from the solution. Moreover, it may not converge, even when the objective function is strongly convex. 

In this library, the stabilization technique proposed in the paper is implemented. It consists in bounding the distance between each pair of successive iterates, which often allows for decreasing the number of BB iterations. When the BB method does not converge, this method makes it convergent. Under suitable assumptions, the method converges globally, despite the fact that no line search is involved, and only gradient values are used. Since the number of stabilization steps is finite, the stabilized version inherits the fast local convergence of the BB method. This technique often allows the BB method to solve problems in a fewer iterations, or even to solve problems where the latter fails.

# Usage

This simple examples shows its most basic use.

````python
>> from scipy.optimize import rosen, rosen_der
>> from stabbb import stab_BB
>> x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
>> bestX, _ = stab_BB(x0, rosen, rosen_der)
>> bestX
array([ 1.,  1.,  1.,  1.,  1.])
````

Please refer to the docstring for a longer explanation.