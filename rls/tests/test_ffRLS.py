#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")


import numpy as np
import pytest
from unittest.mock import patch
from submodules.gcnutils.knavigation import kArrayNav
from submodules.gcnutils.kltisystems import k1OrderLTIsysMimoDiscrete
from rls import kForgettingFactorRLS
from rls.kForgettingFactorRLS import fn_example_scalar, fn_example_vector

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_ffRLS:

    def test_example_scalar(self):
        fn_example_scalar(block=False)

    def test_example_scalar_shall_not_call_update_vector(self):
        with patch.object(kForgettingFactorRLS, "update_vector") as mock:
            mock.return_value = np.nan
            fn_example_scalar(block=False)
            assert mock.call_count == 0

    def test_positional_error(self):
        with pytest.raises(AttributeError):
            rls = kForgettingFactorRLS([0,0,0,0],1.1)

    def test_dimensions_wrong(self):
        with pytest.raises(AssertionError):
            rls = kForgettingFactorRLS(
                    [0,0,0,0],
                    np.asarray([[1,1],[2,2,]]),
                    1.1
            )

    def test_dimensions_right(self):
        rls = kForgettingFactorRLS(
                [0,0,0,0],
                np.eye(4),
                1.1
        )

    def test_example_vector(self):
        fn_example_vector(block=False)

    def test_example_vector_shall_not_call_update_scalar(self):
        with patch.object(kForgettingFactorRLS, "update_scalar") as mock:
            mock.return_value = np.nan
            fn_example_vector(block=False)
            assert mock.call_count == 0

    def test_outputs_of_examples(self):
        tta_scalar = fn_example_scalar(block=False)
        tta_vector = fn_example_vector(block=False)

        for i in range(4):
            assert abs(tta_scalar[i] - tta_vector[i]) < 1e-9


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
