import pytest
from pycovis.postprocess import runtime


def test_version():
    with runtime.Runtime() as pp:
        matlab_version = pp.matlab_version()
        covis_version  = pp.covis_version()
