#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#
# \author    Luciano Augusto Kruk
# \website   www.kruk.eng.br
# \date      2026.0
#
# \description: This file contains an implementation of the 
#               "Error Filtering Online Learning" (EFOL) algorithm.
#
# \license: Please feel free to use and modify this, but keep this header as
#           part of yours. Thanks.
#
# \source:  www.github.com/conconga/estimators
#
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

from submodules.gcnutils.knavigation import kArrayNav
from submodules.gcnutils.kltisystems import (
        k1OrderLTIsysMimoDiscrete,
    )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kEfol:
    """
    Error Filtering Online Learning.
    """

    def __init__(self,
                 filterpole  = -1,
                 dim_alpha   = None,
                 theta0      = None,
                 Ts          = None,
                 Gamma_theta = None,
                 Gamma_error = None,
     ):


        self._ensure_not_a_None(theta0, "theta0")
        self.theta = kArrayNav(theta0, hvector=False)

        self._ensure_not_a_None(Ts, "Ts")
        self.Ts    = Ts

        self._ensure_not_a_None(Gamma_theta, "Gamma_theta")
        self.Gamma_theta = Gamma_theta

        self._ensure_not_a_None(dim_alpha, "dim_alpha")

        if isinstance(filterpole, (int, float, np.float64)):
            filterpole = [float(filterpole)] * dim_alpha

        if isinstance(filterpole, (list, tuple)):
            for p in filterpole:
                if p >= 0:
                    raise ValueError(f'for stability, filterpole shall be < 0')
        else:
            raise ValueError(f'not prepared for type = "{str(type(filterpole))}"')

        if Gamma_error is None:
            self.Gamma_error = kArrayNav(np.eye(dim_alpha))
        else:
            self.Gamma_error = Gamma_error

        # filter L/(s+L) for alpha and beta:
        self.filter = k1OrderLTIsysMimoDiscrete(filterpole, Ts, [0 for i in range(dim_alpha)])


    def _ensure_not_a_None(self, val, val_name):
        if val is None:
            raise ValueError(f'{val_name} shall be defined')

    def _return_a_scalar_if_is_scalar(self, val):
        if isinstance(val, (int, float)):
            return val
        else:
            return val if val.size > 1 else val.reshape(-1)[0]


    def update(self, alpha, beta, hessian):

        if isinstance(alpha, (int, float, list, tuple)):
            alpha = kArrayNav(alpha, hvector=False)

        if isinstance(beta, (int, float, list, tuple)):
            beta = kArrayNav(beta, hvector=False)

        if isinstance(hessian, (list, tuple, np.ndarray)):
            hessian = kArrayNav(hessian)

        assert not any(np.isnan(alpha))
        assert not any(np.isnan(beta))
        assert not any(np.isnan(hessian).reshape(-1))

        # filtered error:
        self.e = kArrayNav(self.filter.update(beta - alpha), hvector=False)
        assert not any(np.isnan(self.e))

        # derivative:
        if hasattr(hessian, "T"):
            ddt = - self.Gamma_theta * hessian.T * self.Gamma_error * self.e
        else:
            # scalar:
            ddt = - self.Gamma_theta * hessian * self.Gamma_error * self.e

        # integration step:
        self.theta += self.Ts * ddt
        assert not any(np.isnan(self.theta))

        return self._return_a_scalar_if_is_scalar(self.theta)

    def get_filtered_error(self):
        return self._return_a_scalar_if_is_scalar(self.e)


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#

