import numpy as np


class PowerLawDistribution:
    """ Based off code from Adam Ginsburg (https://github.com/keflavich/imf)
    """

    def __init__(self, slope, lim):
        self.xmin = lim[0]
        self.xmax = lim[1]
        self.slope = slope

    def pdf(self, x):
        """ Probability distribution function
        """
        if not isinstance(x, np.ndarray):
            x = np.asarray(x)

        if self.slope == -1:
            d = (x**self.slope / np.log(self.xmax / self.xmin)) \
                * (x >= self.xmin) * (x <= self.xmax)
        else:
            d = x**self.slope * (self.slope + 1) \
                / (self.xmax**(self.slope + 1) - self.xmin**(self.slope + 1)) \
                * (x >= self.xmin) * (x <= self.xmax)
        return d

    def ppf(self, x):
        """ Probability point function
        """
        if not isinstance(x, np.ndarray):
            x = np.asarray(x)

        if self.slope == -1:
            p = np.exp(x * np.log(self.xmax)) * self.xmin
        else:
            p = (x * (self.xmax**(self.slope + 1) - self.xmin**(self.slope + 1))
                 + self.xmin**(self.slope + 1))**(1. / (self.slope + 1))
        p[(x < 0) | (x > 1)] = np.nan
        return p

    def sample(self, size):
        """ Sample the distribution function.
        """
        x = np.random.uniform(low=0.0, high=1.0, size=size)
        return self.ppf(x)


class SplitPowerLawDistribution:
    """ Based off code from Adam Ginsburg (https://github.com/keflavich/imf)
    """

    def __init__(self, slopes, limits):
        if type(slopes) != list:
            raise TypeError("slopes must be list")
        if type(limits) != list:
            raise TypeError("limits must be list")
        if len(slopes) != len(limits)-1:
            raise ValueError(
                "Number of slopes not consitent with limits (len(slopes) != len(limits)-1)")

        self.limits = limits
        self.slopes = slopes
        self.xmin = limits[0]
        self.xmax = limits[-1]
        self.powerlaws = []
        for slope, lower_lim, upper_lim in zip(self.slopes, self.limits[:-1], self.limits[1:]):
            self.powerlaws.append(PowerLawDistribution(
                slope, [lower_lim, upper_lim]))

        for ipow, limit in enumerate(self.limits[:-1]):
            if ipow == 0:
                weights = [1.0]
            else:
                current_weight = self.powerlaws[ipow].pdf(
                    limit) / self.powerlaws[ipow-1].pdf(limit)
                weights.append(weights[-1] / current_weight)
            self.weights = np.array(weights) / np.sum(weights)

    def sample(self, size):
        """ Sample the distribution function.
        """
        x = np.random.multinomial(size, self.weights)
        random_value = []

        for xi, powerlaw in zip(x, self.powerlaws):
            if xi > 0:
                random_value.append(powerlaw.sample(xi))

        return np.concatenate(random_value)


def resample_stellar_mass(mform, mmin=0.51, mmax=100, nimf=400, aimf=-2.3):
    """ Resample the stellar masses to create fully sampled IMF. Default parameters correspond to 
        those used in the star-by-star edge simulation.
        In the simulation stellar masses are sampled in mass bins with some resolution (400 bins by
        default). The mass resultion of the IMF is therefore ~0.25 Msol. 
        resambple_imf()
        Given the mass of individual stars in the simulation, this function returns an array of new
        masses for where each stellar mass has been resampled within the bounds of its mass bin. The
        indexing of the returned list correspond to the indexing of the supplied mass array. 
    """
    # Bins of the original initial mass function
    dm = (mmax-mmin)/(nimf-1)
    mbin = np.arange(mmin, mmax+dm, dm)

    mass = np.zeros(len(mform))
    for lower_mbound, upper_mbound in zip(mbin[:-1], mbin[1:]):
        _imf = PowerLawDistribution(
            slope=aimf, lim=[lower_mbound, upper_mbound])
        _mask = (mform >= lower_mbound) & (mform < upper_mbound)
        mass[_mask] = _imf.sample(mform[_mask].size)

    return mass


def sample_stellar_mass_from_spop(number_high_mass, mass_limit=0.51, slope=[-1.3, -2.3, -2.3], limits=[0.08, 0.5, 0.51, 100]):
    """ Sample an IMF from total mass of stellar population to derive stellar properties with same 
        functions as for star-by-star. Default parameters correspond to those used in simulations.
        Provide total number of resolved stars and parameters for broken powerlaw including the individual
        stars as separate bin for normaliation (default parameters).
    """
    _imf = SplitPowerLawDistribution(slope, limits)
    ntot = number_high_mass
    mass = _imf.sample(int(ntot))
    return mass[mass < mass_limit]
