#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")


import numpy as np
import pytest
import matplotlib.pyplot as plt
from submodules.gcnutils.knavigation import kArrayNav
from submodules.gcnutils.kltisystems import k1OrderLTIsysMimoDiscrete
from efol import kEfol


class TestClass_kEfol:
    is_show_with_block = False

    def test_constructor_requires_theta0_Ts_Gamma_and_dim_alpha(self):
        # missing theta0
        with pytest.raises(ValueError):
            kEfol(filterpole=-1, dim_alpha=1, theta0=None, Ts=0.1, Gamma_theta=1.0)

        # missing Ts
        with pytest.raises(ValueError):
            kEfol(filterpole=-1, dim_alpha=1, theta0=[0.0], Ts=None, Gamma_theta=1.0)

        # missing Gamma_theta
        with pytest.raises(ValueError):
            kEfol(filterpole=-1, dim_alpha=1, theta0=[0.0], Ts=0.1, Gamma_theta=None)

        # missing dim_alpha
        with pytest.raises(ValueError):
            kEfol(filterpole=-1, dim_alpha=None, theta0=[0.0], Ts=0.1, Gamma_theta=1.0)


    def test_constructor_rejects_non_negative_filterpole(self):

        # positive pole or zero should raise
        with pytest.raises(ValueError):
            kEfol(filterpole=0.0, dim_alpha=1, theta0=[0.0], Ts=0.1, Gamma_theta=1.0)

        with pytest.raises(ValueError):
            kEfol(filterpole=0.5, dim_alpha=2, theta0=[0.0, 0.0], Ts=0.1, Gamma_theta=1.0)

        # valid negative pole (single) works
        kev = kEfol(filterpole=-2.0, dim_alpha=2, theta0=[0.0, 0.0], Ts=0.1, Gamma_theta=1.0)
        assert hasattr(kev, "filter")


    def test_update_with_scalar_list_and_ndarray_and_kArrayNav(self):
        kev = kEfol(filterpole=-1.0, dim_alpha=2, theta0=[0.0, 0.0], Ts=0.1, Gamma_theta=np.eye(2))

        # inputs as Python lists
        alpha = [0.1, 0.2]
        beta = [0.3, 0.1]

        # simple hessian as identity (or 1D array shaped appropriately)
        hessian = np.eye(2)
        theta_out = kev.update(alpha, beta, hessian)

        # theta should be an object wrapper (kArrayNav) or behave like an
        # array; check not NaN and shape matches
        assert not np.any(np.isnan(np.array(theta_out).reshape(-1)))
        assert np.array(theta_out).size == 2

        # inputs as numpy arrays
        alpha2 = np.array([0.2, 0.2])
        beta2 = np.array([0.1, 0.4])
        theta_out2 = kev.update(alpha2, beta2, hessian)
        assert not np.any(np.isnan(np.array(theta_out2).reshape(-1)))

        # inputs as kArrayNav (assume it accepts a list)
        alpha3 = kArrayNav([0.0, 0.0])
        beta3 = kArrayNav([0.0, 0.1])
        theta_out3 = kev.update(alpha3, beta3, hessian)
        assert not np.any(np.isnan(np.array(theta_out3).reshape(-1)))


    def test_update_raises_on_nan_inputs(self):
        kev = kEfol(filterpole=-1.0, dim_alpha=2, theta0=[0.0, 0.0], Ts=0.1, Gamma_theta=np.eye(2))
        alpha = [np.nan, 0.0]
        beta = [0.1, 0.2]
        hessian = np.eye(2)

        with pytest.raises(AssertionError):
            kev.update(alpha, beta, hessian)

        # NaN in hessian
        alpha = [0.0, 0.0]
        hessian_bad = np.array([[np.nan, 0.0], [0.0, 1.0]])
        with pytest.raises(AssertionError):
            kev.update(alpha, beta, hessian_bad)

    def test_deadzone(self):

        def fn_deadzone(err):
            if abs(err) < 1:
                mask = [0]
            else:
                mask = [1]
            return mask

        # simulation configuration:
        Tmax = 2      # [s]
        Ts   = 1./100 # [Hz]
        T    = np.arange(0,Tmax, Ts)

        # kEfol configuration:
        theta = -1

        efol = kEfol(
                filterpole  = -15,
                dim_alpha   = 1,
                theta0      = theta,
                Ts          = Ts,
                Gamma_theta = 0.5,
                fn_deadzone = fn_deadzone,
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

            if abs(t-0.75) <= Ts:
                assert  0.5 <= abs(beta-alpha) <= 1.0

            if abs(t-1.75) <= Ts:
                assert  0.4 <= abs(beta-alpha) <= 0.5


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
        #.............................................#

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
