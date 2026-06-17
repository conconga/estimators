  #>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

#==================================#
##WWww=--  import section: --=wwWW##

import numpy as np
from   numpy import sin, cos, pi
import matplotlib.pyplot as plt
from efol import kEfol

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class TestClass_Example_2_kEfol:
    is_show_with_block = False

    """
    alfa(t) = 3.3 * sin(2*pi*3*t)
    beta(t) = T0  * cos(2*pi*3*t + T1)   <= approximator

    The goal is to estimate theta = [T0, T1]
    """

    def test_example(self):

        # simulation configuration:
        Tmax = 4      # [s]
        Ts   = 1./100 # [Hz]
        T    = np.arange(0,Tmax, Ts)

        # kEfol configuration:
        theta = np.asarray([1,0]).reshape(2,1)

        efol = kEfol(
                filterpole  = -40,
                dim_alpha   = 1,
                theta0      = theta,
                Ts          = Ts,
                Gamma_theta = np.diag((30, 0.5))
        )

        # setup:
        fbeta = lambda theta,t: theta[0][0] * cos((2*pi*3*t) + theta[1][0])

        # hessian:
        H     = lambda theta,t: [
                            cos((2*pi*3*t) + theta[1][0]),
                            -theta[0][0]*sin((2*pi*3*t) + theta[1][0]),
        ]

        lst_log = list()
        for t in T:

            # known or measured:
            alfa = 3.3 * sin(2*pi*3*t)        # RHS
            alfa += 0.3*np.random.randn()

            # modeled:
            beta = fbeta(theta,t)             # LHS

            # update:
            theta = efol.update(alfa, beta, H(theta,t))
            e     = efol.get_filtered_error()

            # logging:
            lst_log.append((e, alfa, beta, theta.copy()))

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
        ax[0].legend(('alfa', 'beta'))

        M = np.asarray([i[0] for i in lst_log])
        ax[1].plot(T, M)
        ax[1].grid(True)
        ax[1].set_ylabel("filtered error")

        M = np.asarray([i[3].reshape(-1) for i in lst_log])
        ax[2].plot(T, M)
        ax[2].grid(True)
        ax[2].set_ylabel("theta")

        # do not change this order:
        # (this makes matplotlib, Qt and wayland to cooperate!)
        f.canvas.flush_events()
        f.canvas.draw() 

        #.............................................#
        plt.show(block=self.is_show_with_block)
        plt.savefig("example_efol_2.svg", transparent=True)
        #.............................................#

    #====================================#
