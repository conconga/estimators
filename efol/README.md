- [<span class="toc-section-number">1</span>
  Introduction](#introduction)
  - [<span class="toc-section-number">1.1</span> Notation](#notation)
  - [<span class="toc-section-number">1.2</span> In a
    Nutshell](#in-a-nutshell)
  - [<span class="toc-section-number">1.3</span> EFOL](#efol)
    - [<span class="toc-section-number">1.3.1</span> Goal](#goal)
    - [<span class="toc-section-number">1.3.2</span> Error
      Function](#error-function)
    - [<span class="toc-section-number">1.3.3</span> Filtering with
      LTI](#filtering-with-lti)
    - [<span class="toc-section-number">1.3.4</span> Expansion of
      \$(t)\$](#expansion-of-t)
    - [<span class="toc-section-number">1.3.5</span> Control
      Law](#control-law)

# Introduction

The Error Filtering Online Learning (EFOL) was conceived as a scheme to cope
with unkowns not modelled in dynamic systems when the intention is to control
such systems. For example, a dynamic system might have a drift in its behavior
due to temperature, dirt, light, which not always can be modelled. Aiming at
robustness, the alternative in these cases is to enable an adaptive scheme to
compensate the effect of the missing part of the model at the controller.

To my best knowledge, the first time EFOL was published dates 2010, or at least
the first lines of it, by Mr. Marios Polycarpou in "Stable Adaptive Neural
Control Scheme for Nonlinear Systems", but the best reference for it is the
book "Adaptive Approximation Based Control" from Mr. Farrell and Mr.
Polycarpou.

## Notation

- Underlined variables are vectors or matrices.
- A hat will be on the top of approximations or approximators. For example,
    $\hat{f}\left(\hat{\theta}\right)$ has an approximator function $\hat{f}$
    characterized by the approximation parameter $\hat{\theta}$.
- A superscripted triangle points to the best choice of a parameter:
    $\theta^{\triangle}$.
- The vector $\underline{z}(t)$ collects all states and internal variables of
    the system model at instant $t$.

## In a Nutshell

In case you are not familiar with the concept of EFOL, consider a model
represented by a equality, where the LHS of the equality collects what you can
"see" and measure from a system. And the RHS collects what is not known or any
characteristics of the system that changes with time. For example,
$\alpha(t)=0$ represents a system fully described, while $\alpha(t)=\beta(t)$
means a part of the system could not be fully modeled and this uncertainty is
hidden in $\beta(t)$. The trick now is to approximate $\beta(t)$, and adapt its
internal parameters online to have always the best match between $\alpha(t)$
and $\beta(t)$.

## EFOL

### Goal

Let a function $\underline{\hat{f}}$ be the approximation function. It can be
either a linear, or a nonlinear (eg. as a multilayered neural network)
approximator.

Let:

$$
\underline{\alpha}(t) = \underline{\hat{f}}(\underline{z}(t);\underline{\theta}^{\triangle})
                            + \underline{\delta}(t)
$$

with $\underline{\delta}(t)$ being the part of $\underline{\alpha}(t)$ not
covered by $\underline{\hat{f}}$ despite the selection of the best parameters
$\underline{\theta}^{\triangle}$. And let:

$$ \underline{\beta} = \underline{\hat{f}} (\underline{z}(t); \underline{\hat{\theta}}) $$

the estimator dependent on the parameters $\underline{\hat{\theta}}$. As an
online learning process, $\underline{\hat{\theta}}$ will be estimated on the
fly in a way that the error

$$ \underline{u}(t) = \underline{\beta}(t) - \underline{\alpha}(t) $$

continously decreases.

### Error Function

Based on the definition of $\underline{\alpha}(t)$, we shall rewrite $\underline{u}(t)$ as:

$$
\underline{u}(t) = 
    \underline{\hat{f}} (\underline{z}(t); \underline{\hat{\theta}}) -
    \underline{\hat{f}}(\underline{z}(t);\underline{\theta}^{\triangle}) - \underline{\delta}(t)
$$

but we will assign the task to fix the mismatch between $\underline{\alpha}(t)$
and the respective $\underline{\hat{f}}(\cdot)$ to the approximator in
$\underline{\beta}(t)$. Hence:

$$
\underline{u}(t) = 
    \underline{\hat{f}} (\underline{z}(t); \underline{\hat{\theta}}) -
    \underline{\hat{f}}(\underline{z}(t);\underline{\theta}^{\triangle})
$$

### Filtering with LTI

We need a solution that provides both:

- stable and continuous tracking of a signal, and
- a closed algebraic equation for the first-derivative of the tracked signal.

and we can achieve it with Linear Time-Invariant filters.

Using the Laplace notation for LTI filters, let:

$$ \underline{e}(t) = \mathbb{W} [s] ( \underline{u}(t) ) $$

be the signal tracking $\underline{u}(t)$ at the output of the filter $\mathbb{W}[s]$.
The selection of the filter is not important. Now we just need the first derivative, and
therefore we pick this filter:

$$ \mathbb{W}[s] = \frac{\lambda}{s+\lambda} $$

with $\lambda \in \mathbb{R}^{+}$. The differential equation of the output of the filter is:

$$ \underline{\dot{e}}(t) = -\lambda \underline{e}(t) + \lambda \underline{u}(t) $$

### Expansion of $\underline{u}(t)$

Here we present an alternative to have $\underline{u}(t)$ as a function of $\underline{\tilde{\theta}}$:

$$
\begin{aligned}
    \underline{u}(t) & = 
        \underline{\hat{f}} (\underline{z}(t); \underline{\hat{\theta}}) -
        \underline{\hat{f}}(\underline{z}(t);\underline{\theta}^{\triangle}) \\
        & = 
            \underline{\hat{f}} (\underline{z}(t); \underline{\hat{\theta}}) -
            \underline{\hat{f}}(\underline{z}(t);\underline{\hat{\theta}} - \underline{\tilde{\theta}}) \\
        & =
            \underline{\hat{f}} (\underline{z}(t); \underline{\hat{\theta}}) -
            \underline{\hat{f}}(\underline{z}(t);\underline{\hat{\theta}}) +
            \frac{\partial \underline{\hat{f}}}{\partial \underline{\theta}} (\underline{z}(t); \underline{\hat{\theta}}) \cdot
            \underline{\tilde{\theta}}  +  \mathfrak{F} \\
        & \simeq \underline{H}_ {\theta} \cdot \underline{\tilde{\theta}} 
\end{aligned}
$$

### Control Law

The next step is to find a control law for $\underline{\hat{\theta}}$. This is
achieved using Lyapunov to force the absolute value of the filtered error
$\underline{e}(t)$ and the parameters errors $\underline{\tilde{\theta}}$
continuously to decrease. We select the following cost function:

$$ V(t) = \frac{1}{2\lambda} \underline{e}^T \underline{\Gamma}_ e \underline{e} +
            \frac{1}{2} \underline{\tilde{\theta}}^T \underline{\Gamma}_ \theta^{-1} \underline{\tilde{\theta}} $$

where $\underline{\Gamma}_ e$ and $\underline{\Gamma}_ \theta$ are positive semidefinite matrices. After differentiating
the cost function, one gets the control law for $\underline{\hat{\theta}}$:

$$ \underline{\dot{\hat{\theta}}} = - \underline{\Gamma}_ \theta \underline{H}_ {\theta}^T \underline{\Gamma}_ e \underline{e} $$

This equation is implemented in discrete form in the package `kEfol`.



- example

- theory

- how to use it

