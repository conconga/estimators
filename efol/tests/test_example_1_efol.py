  #>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

#==================================#
##WWww=--  import section: --=wwWW##

import numpy as np
import matplotlib.pyplot as plt
from efol import kEfol

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class TestClass_Example_1_kEfol:
    is_show_with_block = False

    """
    Simple test:

    alpha(t) = 16. if t < 1.0 else 9.
    beta(t)  = (2 + theta)^2 <= approximator

    The goal is to estimate theta.
    """

    def test_example(self):

        # simulation configuration:
        Tmax = 2      # [s]
        Ts   = 1./100 # [Hz]
        T    = np.arange(0,Tmax, Ts)

        # kEfol configuration:
        theta = -1

        efol = kEfol(
                filterpole  = -15,
                dim_error   = 1,
                theta0      = theta,
                Ts          = Ts,
                Gamma_theta = 0.5,
        )

        # setup:
        fbeta = lambda theta: (2.+theta)**2 # <= approximator
        H     = lambda theta: 2.*(2+theta)  # <= hessian

        lst_log = list()
        for t in T:

            # known or measured:
            alpha = 16. if t < 1.0 else 9.   # RHS

            # modeled:
            beta = fbeta(theta)             # LHS

            # update:
            theta = efol.update(alpha, beta, H(theta))
            e     = efol.get_filtered_error()

            # logging:
            lst_log.append((e, alpha, beta, theta))

        #.............................................#
        fig = 0

        #---- new figure:
        fig = fig + 1

        f = plt.figure(fig).clf()
        f, ax = plt.subplots(3,1,num=fig)
        ax = ax.reshape(-1)

        M = np.asarray([i[1:3] for i in lst_log])
        ax[0].plot(T, M)
        ax[0].grid(True)
        ax[0].legend(('alpha', 'beta'))

        M = np.asarray([i[0] for i in lst_log])
        ax[1].plot(T, M)
        ax[1].grid(True)
        ax[1].set_ylabel("filtered error")

        M = np.asarray([i[3] for i in lst_log])
        ax[2].plot(T, M)
        ax[2].grid(True)
        ax[2].set_ylabel("theta")

        # do not change this order:
        # (this makes matplotlib, Qt and wayland to cooperate!)
        f.canvas.flush_events()
        f.canvas.draw() 

        #.............................................#
        plt.show(block=self.is_show_with_block)
        plt.savefig("example_efol_1.svg", transparent=True)
        #.............................................#

    #====================================#
