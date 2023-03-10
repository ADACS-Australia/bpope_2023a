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


@pytest.fixture
def model():
    class Model:
        def __init__(self):
            # Instantiate a star with a dipole map
            A = starry.Primary(starry.Map(ydeg=1), prot=0.0)
            amp_true = 0.75
            y_true = np.array([1, 0.1, 0.2, 0.3])
            inc_true = 60
            A.map.amp = amp_true
            A.map[1, :] = y_true[1:]
            A.map.inc = inc_true

            # Instantiate two transiting planets with different longitudes of
            # ascending node. This ensures there's no null space!
            b = starry.Secondary(
                starry.Map(amp=0), porb=1.0, r=0.1, t0=-0.05, Omega=30.0
            )
            c = starry.Secondary(
                starry.Map(amp=0), porb=1.0, r=0.1, t0=0.05, Omega=-30.0
            )
            sys = starry.System(A, b, c)

            # Generate a synthetic light curve with just a little noise
            t = np.linspace(-0.1, 0.1, 100)
            flux = sys.flux(t)
            sigma = 1e-5
            np.random.seed(1)
            flux += np.random.randn(len(t)) * sigma

            # Store
            self.A = A
            self.b = b
            self.c = c
            self.sys = sys
            self.t = t
            self.flux = flux
            self.sigma = sigma
            self.amp_true = amp_true
            self.y_true = y_true
            self.inc_true = inc_true

    return Model()


@pytest.mark.parametrize("L,C", solve_inputs)
def test_solve(L, C, model):
    # Place a generous prior on the map coefficients
    if L == "scalar":
        model.A.map.set_prior(L=1)
    elif L == "vector":
        model.A.map.set_prior(L=np.ones(model.A.map.Ny))
    elif L == "matrix":
        model.A.map.set_prior(L=np.eye(model.A.map.Ny))
    elif L == "cholesky":
        model.A.map.set_prior(cho_L=np.eye(model.A.map.Ny))

    # Provide the dataset
    if C == "scalar":
        model.sys.set_data(model.flux, C=model.sigma ** 2)
    elif C == "vector":
        model.sys.set_data(
            model.flux, C=np.ones(len(model.t)) * model.sigma ** 2
        )
    elif C == "matrix":
        model.sys.set_data(
            model.flux, C=np.eye(len(model.t)) * model.sigma ** 2
        )
    elif C == "cholesky":
        model.sys.set_data(
            model.flux, cho_C=np.eye(len(model.t)) * model.sigma
        )

    # Solve the linear problem
    model.A.map.inc = model.inc_true
    mu, cho_cov = model.sys.solve(t=model.t)

    # Ensure the likelihood of the true value is close to that of
    # the MAP solution
    mu = mu.eval()
    cho_cov = cho_cov.eval()
    cov = cho_cov.dot(cho_cov.T)
    LnL0 = multivariate_normal.logpdf(mu, mean=mu, cov=cov)
    LnL = multivariate_normal.logpdf(
        model.amp_true * model.y_true, mean=mu, cov=cov
    )
    assert LnL0 - LnL < 5.00

    # Check that we can draw from the posterior
    model.sys.draw()


@pytest.mark.parametrize("L,C,woodbury", lnlike_inputs)
def test_lnlike(L, C, woodbury, model):
    # Place a generous prior on the map coefficients
    if L == "scalar":
        model.A.map.set_prior(L=1)
    elif L == "vector":
        model.A.map.set_prior(L=np.ones(model.A.map.Ny))
    elif L == "matrix":
        model.A.map.set_prior(L=np.eye(model.A.map.Ny))
    elif L == "cholesky":
        model.A.map.set_prior(cho_L=np.eye(model.A.map.Ny))

    # Provide the dataset
    if C == "scalar":
        model.sys.set_data(model.flux, C=model.sigma ** 2)
    elif C == "vector":
        model.sys.set_data(
            model.flux, C=np.ones(len(model.t)) * model.sigma ** 2
        )
    elif C == "matrix":
        model.sys.set_data(
            model.flux, C=np.eye(len(model.t)) * model.sigma ** 2
        )
    elif C == "cholesky":
        model.sys.set_data(
            model.flux, cho_C=np.eye(len(model.t)) * model.sigma
        )

    # Compute the marginal log likelihood for different secondari radii
    rs = [0.05, 0.075, 0.1, 0.125, 0.15]
    ll = np.zeros_like(rs)
    for i, r in enumerate(rs):
        model.b.r = r
        ll[i] = model.sys.lnlike(t=model.t, woodbury=woodbury).eval()

    # Verify that we get the correct radius
    assert rs[np.argmax(ll)] == 0.1
    assert np.allclose(ll[np.argmax(ll)], 981.9091)  # benchmarked
