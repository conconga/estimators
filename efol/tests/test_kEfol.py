#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")


import numpy as np
import pytest
from submodules.gcnutils.knavigation import kArrayNav
from submodules.gcnutils.kltisystems import k1OrderLTIsysMimoDiscrete
from efol import kEfol


class TestClass_kEfol:

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

