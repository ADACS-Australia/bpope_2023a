# -*- coding: utf-8 -*-
"""
System linear solve tests.

"""
import starry
import numpy as np
from scipy.linalg import cho_solve
from scipy.stats import multivariate_normal
import pytest
import itertools

# Parameter combinations we'll test
# This test is slow; since it's redundant with
# `test_solve_lazy`, let's only do a subset here
vals = ["scalar"]  # ["vector", "matrix", "cholesky"]
woodbury = [False, True]
solve_inputs = itertools.product(vals, vals)
lnlike_inputs = itertools.product(vals, vals, woodbury)

# Instantiate a star with a dipole map
A = starry.Primary(starry.Map(ydeg=1), prot=0.0)
y_true = [0.1, 0.2, 0.3]
inc_true = 60
A.map[1, :] = y_true
A.map.inc = inc_true

# Instantiate two transiting planets with different longitudes of
# ascending node. This ensures there's no null space!
b = starry.Secondary(starry.Map(), porb=1.0, r=0.1, t0=-0.05, Omega=30.0)
c = starry.Secondary(starry.Map(), porb=1.0, r=0.1, t0=0.05, Omega=-30.0)
sys = starry.System(A, b, c)

# Generate a synthetic light curve with just a little noise
t = np.linspace(-0.1, 0.1, 100)
flux = sys.flux(t)
sigma = 1e-5
np.random.seed(1)
flux += np.random.randn(len(t)) * sigma


@pytest.mark.parametrize("L,C", solve_inputs)
def test_solve(L, C):
    # Place a generous prior on the map coefficients
    if L == "scalar":
        A.map.set_prior(L=1)
    elif L == "vector":
        A.map.set_prior(L=np.ones(A.map.Ny - 1))
    elif L == "matrix":
        A.map.set_prior(L=np.eye(A.map.Ny - 1))
    elif L == "cholesky":
        A.map.set_prior(cho_L=np.eye(A.map.Ny - 1))

    # Provide the dataset
    if C == "scalar":
        sys.set_data(flux, C=sigma ** 2)
    elif C == "vector":
        sys.set_data(flux, C=np.ones(len(t)) * sigma ** 2)
    elif C == "matrix":
        sys.set_data(flux, C=np.eye(len(t)) * sigma ** 2)
    elif C == "cholesky":
        sys.set_data(flux, cho_C=np.eye(len(t)) * sigma)

    # Solve the linear problem
    A.map.inc = inc_true
    mu, cho_cov = sys.solve(t=t)

    # Ensure the likelihood of the true value is close to that of
    # the MAP solution
    mean = mu[0].eval()
    cho_cov = cho_cov[0].eval()
    cov = cho_cov.dot(cho_cov.T)
    LnL0 = multivariate_normal.logpdf(mean, mean=mean, cov=cov)
    LnL = multivariate_normal.logpdf([0.1, 0.2, 0.3], mean=mean, cov=cov)
    assert LnL0 - LnL < 5.00


@pytest.mark.parametrize("L,C,woodbury", lnlike_inputs)
def test_lnlike(L, C, woodbury):
    # Place a generous prior on the map coefficients
    if L == "scalar":
        A.map.set_prior(L=1)
    elif L == "vector":
        A.map.set_prior(L=np.ones(A.map.Ny - 1))
    elif L == "matrix":
        A.map.set_prior(L=np.eye(A.map.Ny - 1))
    elif L == "cholesky":
        A.map.set_prior(cho_L=np.eye(A.map.Ny - 1))

    # Provide the dataset
    if C == "scalar":
        sys.set_data(flux, C=sigma ** 2)
    elif C == "vector":
        sys.set_data(flux, C=np.ones(len(t)) * sigma ** 2)
    elif C == "matrix":
        sys.set_data(flux, C=np.eye(len(t)) * sigma ** 2)
    elif C == "cholesky":
        sys.set_data(flux, cho_C=np.eye(len(t)) * sigma)

    # Compute the marginal log likelihood for different secondari radii
    rs = [0.05, 0.075, 0.1, 0.125, 0.15]
    ll = np.zeros_like(rs)
    for i, r in enumerate(rs):
        b.r = r
        ll[i] = sys.lnlike(t=t, woodbury=woodbury).eval()

    # Verify that we get the correct radius
    assert rs[np.argmax(ll)] == 0.1
    assert np.allclose(ll[np.argmax(ll)], 989.3974)  # benchmarked