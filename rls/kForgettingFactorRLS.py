#~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-
#
# \author    Luciano Augusto Kruk
# \website   www.kruk.eng.br
# \date      2017.0
#
# \description: This file contains an implementation of the forgetting-factor
#               RLS algorithm, with an unit test of first order continuos
#               transfer function.
#
# \license: Please feel free to use and modify this, but keep this header as
#           part of yours. Thanks.
#
# \source:  www.github.com/conconga/estimators
#
#~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-
##WWww=--  import section: --=wwWW##
import numpy               as np

###########################
## Forgetting Factor RLS ##
###########################
class kForgettingFactorRLS:
    """
    Estimation of the vector \\theta when either the output is a scalar value or a vector.

    For y scalar:
        y[k] = h[k]' . theta

    or for y a vector:
        y[k] = H[k] . theta

    where:
        
        theta  :   [N x 1] vector
        h[k]   :   [N x 1] vector
        H[k]   :   [M x N] matrix

        y[k]   :   (float)
    or
        y[k]   :   [M x 1] vector

    """

    def __init__(self, theta_0, P_0=None, lbd=1.):
        """
        Inputs:
            
            theta_0   :    [N x 1] initial guess of the coeficients
            P_0       :    [N x N] initial guess of cov matrix
            lbd       :    (float) forgetting factor (0<lbd<=1)
        """

        self.N       = len(theta_0)
        self.lbd     = lbd 
        self.theta   = np.asarray(theta_0).reshape((self.N,1))

        if P_0 is None:
            self.P = np.eye(self.N)
        else:
            assert P_0.shape == (self.N, self.N), ValueError('matrix P shall be [N x N]')
            self.P = P_0

    def update(self, y, hH):

        if isinstance(y, (int, float)):
            return self.update_scalar(y, hH)

        else:
            return self.update_vector(y, hH)


    def update_scalar(self, y, h):
        """
        y = h^T . theta
        """

        h = h.reshape((self.N,1))

        gamma       = (self.P @ h) / (self.lbd + (h.T @ self.P @ h))
        self.theta  = self.theta + (gamma @ (y - (h.T @ self.theta)))
        self.P      = (self.P - (gamma * (h.T @ self.P))) / self.lbd

        return self.theta

    def update_vector(self, y, H):
        """
        y = H . theta
        """

        M,N = H.shape

        K           = (self.P @ H.T) @ np.linalg.inv( (self.lbd*np.eye(M)) + (H @ self.P @ H.T) )
        self.theta += K @ (y - (H @ self.theta))
        self.P      = (1./self.lbd) * ((np.eye(N) - (K @ H)) @ self.P)

        return self.theta

#################################
def fn_example_scalar(block=False):

    from submodules.gcnutils.kltisystems.k1orderltisyssiso import (
            k1OrderLTIsysSisoDiscrete as CF1ORD_D,
    )

    import matplotlib.pyplot  as plt

    a    = -15.  # pole for the transfer function:
    tmax = 2.0
    Ts   = 5e-3 # sample rate
    T    = [i*Ts for i in range(int(tmax/Ts)+1)] # time vector
    sys  = CF1ORD_D(a, Ts, 1) # unknown system
    rls  = kForgettingFactorRLS(np.random.randn(4), lbd=0.85)

    log_h     = [0]
    log_y     = [sys.y]
    log_theta = [rls.theta.squeeze().tolist()]

    for t in T:
        if (not (t%0.5)):
            h = float(np.random.rand()) # plant input

        y = sys.update(h)
        rls.update(-0.55+y, np.asarray([
            1.,
            h,
            log_h[-1],
            log_y[-1]
        ]))
        
        log_h.append(float(h))
        log_y.append(sys.y)
        log_theta.append(rls.theta.squeeze().tolist())

    T.append(tmax)
    

    #.............................................#
    #---- new figure:
    fig = 0

    #---- new figure:
    fig = fig + 1

    f = plt.figure(fig).clf()
    f, ax = plt.subplots(2,1,num=fig)

    ax[0].plot(T, log_h, T, log_y)
    ax[0].grid(True)

    ax[1].plot(
            T, log_theta, 
            plt.xlim(), (sys.b, sys.b), 'k--', 
            plt.xlim(), (sys.c, sys.c), 'k--',
            plt.xlim(), (-0.55, -0.55), 'k--'
    )

    ax[1].grid(True, alpha=0.6)

    #.............................................#
    # (this makes matplotlib, Qt and wayland to cooperate!)
    for fig in plt.get_fignums():
        plt.figure(fig).canvas.flush_events()
        plt.figure(fig).canvas.draw() 

    plt.show(block=block)
    #.............................................#
    return log_theta[-1]

#~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-

def fn_example_vector(block=False):
    from submodules.gcnutils.kltisystems.k1orderltisysmimo import (
            k1OrderLTIsysMimoDiscrete as MIMO,
    )

    import matplotlib.pyplot  as plt

    a    = -15.  # pole for the transfer function:
    tmax = 2.0
    Ts   = 5e-3 # sample rate
    T    = np.arange(0,tmax,Ts)
    mimo = MIMO(a, Ts, [0,0])
    rls  = kForgettingFactorRLS(np.random.randn(4), lbd=0.85)

    log_in  = list()
    log_out = list()
    log_tta = list()
    inputs  = np.zeros(2)

    for t in T:
        last_inputs = inputs.copy()
        if (not (t%0.5)):
            inputs = np.random.rand(2)

        last_outputs = mimo.y.copy()
        y = mimo.update(inputs).reshape(2,1)

        rls.update(-0.55+y, np.asarray([
            [ 1., inputs[0], last_inputs[0], last_outputs[0]],
            [ 1., inputs[1], last_inputs[1], last_outputs[1]],
        ]))

        log_in.append(inputs)
        log_out.append(mimo.y)
        log_tta.append(rls.theta.squeeze().tolist())

    # figure:
    nb_fig = 10
    plt.figure(nb_fig).clf()
    fig, ax = plt.subplots(2,1,num=nb_fig,sharex=True)

    ax[0].plot(T, np.asarray(log_in), T, np.asarray(log_out))
    ax[0].grid(True, alpha=0.6)

    ax[1].plot(T, np.asarray(log_tta))
    ax[1].grid(True, alpha=0.6)


    #.............................................#
    # (this makes matplotlib, Qt and wayland to cooperate!)
    for fig in plt.get_fignums():
        plt.figure(fig).canvas.flush_events()
        plt.figure(fig).canvas.draw() 

    plt.show(block=block)
    #.............................................#
    return log_tta[-1]

#~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-

if __name__=="__main__":
    fn_example_scalar(block=False)
    fn_example_vector(block=False)

#~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-~-=-
