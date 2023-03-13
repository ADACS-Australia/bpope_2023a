def X(
    self,
    t,
    pri_r,
    pri_m,
    pri_prot,
    pri_t0,
    pri_theta0,
    pri_amp,
    pri_inc,
    pri_obl,
    pri_fproj,
    pri_u,
    pri_f,
    sec_r,
    sec_m,
    sec_prot,
    sec_t0,
    sec_theta0,
    sec_porb,
    sec_ecc,
    sec_w,
    sec_Omega,
    sec_iorb,
    sec_amp,
    sec_inc,
    sec_obl,
    sec_u,
    sec_f,
    sec_sigr,
):
    """Compute the system light curve design matrix."""
    # # Exposure time integration?
    if self.texp != 0.0:
        print("system x: integration")
        pass
        
    # Compute the relative positions of all bodies
    orbit = exoplanet.orbits.KeplerianOrbit(
        period=sec_porb,
        t0=sec_t0,
        incl=sec_iorb,
        ecc=sec_ecc,
        omega=sec_w,
        Omega=sec_Omega,
        m_planet=sec_m,
        m_star=pri_m,
        r_star=pri_r,
    )
    try:
        x, y, z = orbit.get_relative_position(
            t, light_delay=self.light_delay
        )
    except TypeError:
        if self.light_delay:
            logger.warn(
                "This version of `exoplanet` does not model light delays."
            )
        x, y, z = orbit.get_relative_position(t)

    # The shape of `x`, `y`, `z` in `exoplanet` is unreliable!
    if x.ndim == 1:
        x = tt.reshape(x, (-1, 1))
        y = tt.reshape(y, (-1, 1))
        z = tt.reshape(z, (-1, 1))

    # Get all rotational phases
    pri_prot = ifelse(
        tt.eq(pri_prot, 0.0), math.to_tensor(np.inf), pri_prot
    )
    print("pri_prot", pri_prot)
    theta_pri = (2 * np.pi) / pri_prot * (t - pri_t0) + pri_theta0
    sec_prot = tt.switch(
        tt.eq(sec_prot, 0.0), math.to_tensor(np.inf), sec_prot
    )
    print("sec prot: ", sec_prot)
    theta_sec = (2 * np.pi) / tt.shape_padright(sec_prot) * (
        tt.shape_padleft(t) - tt.shape_padright(sec_t0)
    ) + tt.shape_padright(sec_theta0)

    # Compute all the phase curves
    if self._oblate:
        pass
    else:
        print("oblate false")
        phase_pri = pri_amp * self.primary.map.ops.X(
            theta_pri,
            tt.zeros_like(t),
            tt.zeros_like(t),
            tt.zeros_like(t),
            math.to_tensor(0.0),
            pri_inc,
            pri_obl,
            pri_u,
            pri_f,
        )
    if self._reflected:
        pass
    else:
        print("reflected false")
        phase_sec = [
            sec_amp[i]
            * sec.map.ops.X(
                theta_sec[i],
                -x[:, i],
                -y[:, i],
                -z[:, i],
                math.to_tensor(0.0),  # occultor of zero radius
                sec_inc[i],
                sec_obl[i],
                sec_u[i],
                sec_f[i],
            )
            for i, sec in enumerate(self.secondaries)
        ]

    # Compute any occultations
    occ_pri = tt.zeros_like(phase_pri)
    occ_sec = [tt.zeros_like(ps) for ps in phase_sec]

    # Compute the period if we were given a semi-major axis
    sec_porb = tt.switch(
        tt.eq(sec_porb, 0.0),
        (G_grav * (pri_m + sec_m) * sec_porb ** 2 / (4 * np.pi ** 2))
        ** (1.0 / 3),
        sec_porb,
    )

    # Compute transits across the primary
    for i, _ in enumerate(self.secondaries):
        xo = x[:, i] / pri_r
        yo = y[:, i] / pri_r
        zo = z[:, i] / pri_r
        ro = sec_r[i] / pri_r
        b = tt.sqrt(xo ** 2 + yo ** 2)
        b_occ = tt.invert(
            tt.ge(b, 1.0 + ro) | tt.le(zo, 0.0) | tt.eq(ro, 0.0)
        )
        idx = tt.arange(b.shape[0])[b_occ]
        if self._oblate:
            pass
        else:
            occ_pri = tt.set_subtensor(
                occ_pri[idx],
                occ_pri[idx]
                + pri_amp
                * self.primary.map.ops.X(
                    theta_pri[idx],
                    xo[idx],
                    yo[idx],
                    zo[idx],
                    ro,
                    pri_inc,
                    pri_obl,
                    pri_u,
                    pri_f,
                )
                - phase_pri[idx],
            )

    # Compute occultations by the primary
    for i, sec in enumerate(self.secondaries):

        xo = -x[:, i] / sec_r[i]
        yo = -y[:, i] / sec_r[i]
        zo = -z[:, i] / sec_r[i]
        ro = pri_r / sec_r[i]
        b = tt.sqrt(xo ** 2 + yo ** 2)
        b_occ = tt.invert(
            tt.ge(b, 1.0 + ro) | tt.le(zo, 0.0) | tt.eq(ro, 0.0)
        )
        idx = tt.arange(b.shape[0])[b_occ]
        if self._oblate:
            # TODO: Occultations *by* an oblate occultor are not
            # currently supported. The following code ignores any
            # oblateness and instead treats the body as a spherical
            # occultor with radius equal to its equatorial radius.
            pass
        elif self._reflected:
            pass
        else:
            occ_sec[i] = tt.set_subtensor(
                occ_sec[i][idx],
                occ_sec[i][idx]
                + sec_amp[i]
                * sec.map.ops.X(
                    theta_sec[i, idx],
                    xo[idx],
                    yo[idx],
                    zo[idx],
                    ro,
                    sec_inc[i],
                    sec_obl[i],
                    sec_u[i],
                    sec_f[i],
                )
                - phase_sec[i][idx],
            )

    # Compute secondary-secondary occultations
    for i, sec in enumerate(self.secondaries):
        for j, _ in enumerate(self.secondaries):
            if i == j:
                continue
            xo = (-x[:, i] + x[:, j]) / sec_r[i]
            yo = (-y[:, i] + y[:, j]) / sec_r[i]
            zo = (-z[:, i] + z[:, j]) / sec_r[i]
            ro = sec_r[j] / sec_r[i]
            b = tt.sqrt(xo ** 2 + yo ** 2)
            b_occ = tt.invert(
                tt.ge(b, 1.0 + ro) | tt.le(zo, 0.0) | tt.eq(ro, 0.0)
            )
            idx = tt.arange(b.shape[0])[b_occ]
            if self._reflected:
                pass
            else:
                occ_sec[i] = tt.set_subtensor(
                    occ_sec[i][idx],
                    occ_sec[i][idx]
                    + sec_amp[i]
                    * sec.map.ops.X(
                        theta_sec[i, idx],
                        xo[idx],
                        yo[idx],
                        zo[idx],
                        ro,
                        sec_inc[i],
                        sec_obl[i],
                        sec_u[i],
                        sec_f[i],
                    )
                    - phase_sec[i][idx],
                )

    # Concatenate the design matrices
    X_pri = phase_pri + occ_pri
    X_sec = [ps + os for ps, os in zip(phase_sec, occ_sec)]
    X = tt.horizontal_stack(X_pri, *X_sec)

    # Sum and return
    if self.texp == 0.0:
        print("x texp 0")

        return X

    else:
        pass




#### OpsLD
def X(self, theta, xo, yo, zo, ro, inc, obl, u, f):
    """
    Convenience function for integration of limb-darkened maps
    with the ``System`` class. The design matrix for limb-darkened
    maps is just a column vector equal to the total flux, since the
    spherical harmonic coefficient vector is ``[1.0]``.

    """
    print("es call: LD_X")
    flux = self.flux(xo, yo, zo, ro, u)
    X = tt.reshape(flux, (-1, 1))
    return X
