---
layout: post
title: A double pendulum subjected to a follower load
date: 2023-05-07
description: Numerical implementation of the stability analysis of a follower force system.
categories: numerical-implementation mechanics
related_posts: false
---

## Governing equations

The derivation in this section generally follows the method in {% cite  Semler1998 --file bibliography %}.

The system, under zero gravity, consists of two masses, $$2m$$ and $$m$$, interconnected by massless rigid bars of length $$L$$, the lower one of which is subjected to a follower load, $$P$$.
The two bars are constrained by rotational springs of equal stiffness, $$k$$, and rotational dashpots, $$c_1$$ and $$c_2$$.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/posts/Follower_Force_System.png" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Sketch of the double pendulum system in zero gravity with a follower load
</div>

Introducing the dimensionless time $$\tau=t\sqrt{k/(mL^2)}$$ and the following dimensionless parameters
\begin{equation}
\label{eq:dimensionless_parameters}
\gamma_1 = \frac{c_1}{\sqrt{kmL^2}},\quad \gamma_2 = \frac{c_2}{\sqrt{kmL^2}},\quad p=\frac{PL}{k},
\end{equation}
the dimensionless governing equation will be
\begin{equation}
\label{eq:dimensionless_ge}
\boldsymbol{M}\ddot{\boldsymbol{\phi}} + \boldsymbol{C}\dot{\boldsymbol{\phi}} + \boldsymbol{K}\boldsymbol{\phi} = 0
\end{equation}
where

$$
\boldsymbol{M} = 
\begin{bmatrix}
3 & 1\\
1 & 1
\end{bmatrix},\quad
\boldsymbol{C} = 
\begin{bmatrix}
\gamma_1 + \gamma_2 & -\gamma_2\\
-\gamma_2 & \gamma_2
\end{bmatrix},\quad
\boldsymbol{K} = 
\begin{bmatrix}
2 - p & p - 1\\
-1 & 1
\end{bmatrix}
$$

and

$$
\boldsymbol{\phi} = 
\begin{bmatrix}
\phi_1 \\
\phi_2
\end{bmatrix}
$$

Now it is clear that the system is defined by three parameters, the two nondimensional dmaping coefficients, $$\gamma_1$$ and $$\gamma_2$$, and the nondimensional follower force, $$p$$.

Furthermore, defining the following partitioned vector and matrices of order $$4$$

$$
\boldsymbol{\varphi} = 
\begin{bmatrix}
\dot{\boldsymbol{\phi}}\\
\boldsymbol{\phi}
\end{bmatrix},\quad
\boldsymbol{B} = 
\begin{bmatrix}
\boldsymbol{0} & \boldsymbol{M} \\
\boldsymbol{M} & \boldsymbol{C} 
\end{bmatrix},\quad
\boldsymbol{E} = 
\begin{bmatrix}
-\boldsymbol{M} & \boldsymbol{0} \\
\boldsymbol{0} & \boldsymbol{K} 
\end{bmatrix}
$$

the problem reduces into the first-order form
\begin{equation}
\label{eq:first_order_ge}
\boldsymbol{B}\dot{\boldsymbol{\varphi}} + \boldsymbol{E}\boldsymbol{\varphi} = 0
\end{equation}

## Eigenvalue problem

Assuming solutions of the periodic form $$\boldsymbol{\varphi}=\boldsymbol{A}\exp(\lambda t) \equiv \boldsymbol{A}\exp(i\omega t)$$,
where $$\boldsymbol{A}$$ is the modal amplitude and $$i=\sqrt{-1}$$, leads to the eigenvalue problem
\begin{equation}
\label{eq:eigenvalue_problem}
(\lambda\boldsymbol{I} - \boldsymbol{Y}) \boldsymbol{A} = \boldsymbol{0}
\end{equation}
where $$\boldsymbol{Y}=-\boldsymbol{B}^{-1}\boldsymbol{E}$$.
The eigenvalues $$\lambda_i$$ are in complex and $$\lambda_i=i\omega$$; the real part of $$\omega$$  is associated with frequency of oscillation and the imaginary part with damping.

Since $$\boldsymbol{K}$$ can be nonsymmetrical (when $$p\ne 0$$), $$\boldsymbol{Y}$$ may be nonsymmetrical as well, which is a general eigenvalue problem{% cite  Meirovitch1967 --file bibliography %}.
In addition to Eq. $$\eqref{eq:eigenvalue_problem}$$, the adjoint of eigenvalue problem is defined
\begin{equation}
\label{eq:adjoint_eigenvalue_problem}
(\lambda\boldsymbol{I} - \boldsymbol{Y}^T) \tilde{\boldsymbol{A}} = \boldsymbol{0}
\end{equation}
the eigenvalues of which are the same as those of Eq. $$\eqref{eq:eigenvalue_problem}$$, but the eigenvectors, $$\tilde{\boldsymbol{A}}$$, are different.
Introducing the transformation $$\boldsymbol{\varphi}=\boldsymbol{A}\boldsymbol{y}$$, the original problem, Eq. $$\eqref{eq:first_order_ge}$$ can be decoupled into
\begin{equation}
\boldsymbol{y} = \boldsymbol{D}\dot{\boldsymbol{y}}
\end{equation}
where $$\boldsymbol{D}=\boldsymbol{P}^{-1}\boldsymbol{S}$$ is diagonal, $$\boldsymbol{P}=\tilde{\boldsymbol{A}}^T\boldsymbol{B}\boldsymbol{A}$$, and $$\boldsymbol{S}=-\tilde{\boldsymbol{A}}^T\boldsymbol{E}\boldsymbol{A}$$.

Using the model matrix $$\boldsymbol{A}$$ and the coordinate transformation, is is then possible to calculate the solution at any time in the original coordinates, $$\phi_1$$ and $$\phi_2$$.


## Dynamic instability

The critical condition for the follower load is directly given
\begin{equation}
    p_cr = \frac{4\gamma_1^2 + 33\gamma_1\gamma_2 + 4\gamma_2^2}{2(\gamma_1^2 + 7\gamma_1\gamma_2 + 6\gamma_2^2)} + \frac{1}{2}\gamma_1\gamma_2
\end{equation}
See {% cite  Semler1998 --file bibliography %} for more details.

## Numerical implementation

A ``Python`` script is given to implement the analytical solution drived before.

```python
import numpy as np

def critical_p(gamma_1, gamma_2):
    return (4 * gamma_1 * gamma_1 + 33 * gamma_1 * gamma_2 + 4 * gamma_2 * gamma_2) / (2 * (gamma_1 * gamma_1 + 7 * gamma_1 * gamma_2 + 6 * gamma_2 * gamma_2)) + 1 / 2 * gamma_1 * gamma_2

# governing parameters
gamma_1 = 0.10
gamma_2 = 0.02
p = 2.0
p_cr = critical_p(gamma_1, gamma_2)
print("gamma_1={:.2f}, gamma_2={:.2f}, p={:.2f}, critical p={:.2f}".format(gamma_1, gamma_2, p, p_cr))
if p > p_cr:
    print("For given conditions, the system is UNSTEABLE!")
else:
    print("For given conditions, the system is STEABLE!")
p = p_cr

# eigenvalue problem
## define matrices
M = np.matrix([[3, 1], [1, 1]])
C = np.matrix([[gamma_1 + gamma_2, -gamma_2], [-gamma_2, gamma_2]])
K = np.matrix([[2 - p, p - 1], [-1, 1]])

B = np.matrix(np.zeros([4, 4]))
E = np.matrix(np.zeros([4, 4]))
B[:2, 2:] =  M
B[2:, :2] =  M
B[2:, 2:] =  C
E[:2, :2] = -M
E[2:, 2:] =  K

Y = -np.linalg.inv(B) * E

## solve the eigenvalue problem and the adjoint one
u1, A1 = np.linalg.eig(Y)
u2, A2 = np.linalg.eig(Y.T)

## sort the eigenvalues and eigenvectors in the same order
def argsort(a):
    # create a new array of tuples that pairs each element with its index
    indexed_a = [(x, i) for i, x in enumerate(a)]

    # define a custom sorting key that sorts by real part first, and then imaginary part
    def complex_sort_key(x):
        # return (x[0].imag, x[0].real)
        return (x[0].real, x[0].imag)

    # sort the indexed array using the custom sorting key
    sorted_indexed_a = sorted(indexed_a, key=complex_sort_key)

    # extract the sorted indices from the sorted indexed array
    sorted_indices = [x[1] for x in sorted_indexed_a]
    
    return sorted_indices

ind1 = argsort(u1)
u1   = u1[ind1]
A1   = A1[:, ind1]
ind2 = argsort(u2)
u2   = u2[ind2]
A2   = A2[:, ind2]

## decoupled problem
P = A2.T * B * A1
S = A2.T * E * A1
D = np.diag(-np.linalg.inv(P) * S)

# define initial conditions
## For example, the double pendulum is released from static horizontal position
phi_1_0 = np.pi / 6
phi_2_0 = np.pi / 3
phi_1_dt_0 = 0
phi_2_dt_0 = 0
## vector in original coordinates
phi_0 = np.matrix([[phi_1_dt_0], [phi_2_dt_0], [phi_1_0], [phi_2_0]])
## vector in transformed coordinates
y0 = np.linalg.inv(A1) * phi_0

def solution(t):
    y = np.exp(u1 * t) * np.reshape(np.array(y0), (1, 4))
    return np.array((A1 * np.matrix(y).T)).real
```

## Reference

{% bibliography --cited --file bibliography %}