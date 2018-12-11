import numexpr
import flarestack.core.astro
import numpy as np
import scipy.interpolate
import cPickle as Pickle
from flarestack.shared import acceptance_path
from flarestack.core.time_PDFs import TimePDF
from flarestack.utils.make_SoB_splines import load_spline, \
    load_bkg_spatial_spline
from flarestack.core.energy_PDFs import EnergyPDF
from flarestack.utils.create_acceptance_functions import dec_range
from flarestack.utils.dataset_loader import data_loader
from flarestack.utils.make_SoB_splines import create_2d_ratio_spline


class LLH:
    """Base class LLH.
    """
    subclasses = {}

    def __init__(self, season, sources, llh_dict):
        self.season = season
        self.sources = sources

        try:
            time_dict = llh_dict["LLH Time PDF"]
            self.time_pdf = TimePDF.create(time_dict, season)
        except KeyError:
            raise KeyError("No Time PDF specified. Please add an 'llh_time_pdf'"
                           " entry to the llh_dict, and try again. If "
                           "you do not want time dependence in your "
                           "likelihood, please specify a 'Steady' Time PDF.")

        self.bkg_spatial = load_bkg_spatial_spline(self.season)
        self.acceptance_f = self.create_acceptance_function()

    @classmethod
    def register_subclass(cls, llh_name):
        """Adds a new subclass of EnergyPDF, with class name equal to
        "energy_pdf_name".
        """
        def decorator(subclass):
            cls.subclasses[llh_name] = subclass
            return subclass

        return decorator

    @classmethod
    def create(cls, season, sources, llh_dict):
        llh_name = llh_dict["name"]

        if llh_name not in cls.subclasses:
            raise ValueError('Bad LLH name {}'.format(llh_name))

        return cls.subclasses[llh_name](season, sources, llh_dict)

    @classmethod
    def get_parameters(cls, llh_dict):
        llh_name = llh_dict["name"]

        if llh_name not in cls.subclasses:
            raise ValueError('Bad LLH name {}'.format(llh_name))

        return cls.subclasses[llh_name].return_llh_parameters(llh_dict)

    @classmethod
    def get_injected_parameters(cls, mh_dict):
        llh_name = mh_dict["llh_dict"]["name"]

        if llh_name not in cls.subclasses:
            raise ValueError('Bad LLH name {}'.format(llh_name))

        return cls.subclasses[llh_name].return_injected_parameters(mh_dict)

    # ==========================================================================
    # Signal PDF
    # ==========================================================================

    def signal_pdf(self, source, cut_data):
        """Calculates the value of the signal spatial PDF for a given source
        for each event in the coincident data subsample. If there is a Time PDF
        given, also calculates the value of the signal Time PDF for each event.
        Returns either the signal spatial PDF values, or the product of the
        signal spatial and time PDFs.

        :param source: Source to be considered
        :param cut_data: Subset of Dataset with coincident events
        :return: Array of Signal Spacetime PDF values
        """
        space_term = self.signal_spatial(source, cut_data)

        if hasattr(self, "time_pdf"):
            time_term = self.time_pdf.signal_f(
                cut_data[self.season["MJD Time Key"]], source)

            sig_pdf = space_term * time_term

        else:
            sig_pdf = space_term

        return sig_pdf

    @staticmethod
    def signal_spatial(source, cut_data):
        """Calculates the angular distance between the source and the
        coincident dataset. Uses a Gaussian PDF function, centered on the
        source. Returns the value of the Gaussian at the given distances.

        :param source: Single Source
        :param cut_data: Subset of Dataset with coincident events
        :return: Array of Spatial PDF values
        """
        distance = flarestack.core.astro.angular_distance(
            cut_data['ra'], cut_data['dec'], source['ra'], source['dec'])
        space_term = (1. / (2. * np.pi * cut_data['sigma'] ** 2.) *
                      np.exp(-0.5 * (distance / cut_data['sigma']) ** 2.))
        return space_term

    # ==========================================================================
    # Background PDF
    # ==========================================================================

    def background_pdf(self, source, cut_data):
        """Calculates the value of the background spatial PDF for a given
        source for each event in the coincident data subsample. Thus is done
        by calling the self.bkg_spline spline function, which was fitted to
        the Sin(Declination) distribution of the data.

        If there is a signal Time PDF given, then the background time PDF
        is also calculated for each event. This is assumed to be a normalised
        uniform distribution for the season.

        Returns either the background spatial PDF values, or the product of the
        background spatial and time PDFs.

        :param source: Source to be considered
        :param cut_data: Subset of Dataset with coincident events
        :return: Array of Background Spacetime PDF values
        """
        space_term = self.background_spatial(cut_data)

        if hasattr(self, "time_pdf"):
            time_term = self.time_pdf.background_f(
                cut_data[self.season["MJD Time Key"]], source)

            sig_pdf = space_term * time_term
        else:
            sig_pdf = space_term

        return sig_pdf

    def background_spatial(self, cut_data):
        space_term = (1. / (2. * np.pi)) * np.exp(
            self.bkg_spatial(cut_data["sinDec"]))
        # space_term = (1 / (4 * np.pi))
        return space_term

    def acceptance(self, source, params=None):
        """Calculates the detector acceptance for a given source, using the
        1D interpolation of the acceptance as a function of declination based
        on the IC data rate. This is a crude estimation.

        :param source: Source to be considered
        :param params: Parameter array
        :return: Value for the acceptance of the detector, in the given
        season, for the source
        """
        return self.acceptance_f(source["dec"])

    def create_acceptance_function(self):
        pass

    def select_spatially_coincident_data(self, data, sources):
        """Checks each source, and only identifies events in data which are
        both spatially and time-coincident with the source. Spatial
        coincidence is defined as a +/- 5 degree box centered on the  given
        source. Time coincidence is determined by the parameters of the LLH
        Time PDF. Produces a mask for the dataset, which removes all events
        which are not coincident with at least one source.

        :param data: Dataset to be tested
        :param sources: Sources to be tested
        :return: Mask to remove
        """
        veto = np.ones_like(data["ra"], dtype=np.bool)

        for source in sources:

            # Sets half width of spatial box
            width = np.deg2rad(5.)

            # Sets a declination band 5 degrees above and below the source
            min_dec = max(-np.pi / 2., source['dec'] - width)
            max_dec = min(np.pi / 2., source['dec'] + width)

            # Accepts events lying within a 5 degree band of the source
            dec_mask = np.logical_and(np.greater(data["dec"], min_dec),
                                      np.less(data["dec"], max_dec))

            # Sets the minimum value of cos(dec)
            cos_factor = np.amin(np.cos([min_dec, max_dec]))

            # Scales the width of the box in ra, to give a roughly constant
            # area. However, if the width would have to be greater that +/- pi,
            # then sets the area to be exactly 2 pi.
            dPhi = np.amin([2. * np.pi, 2. * width / cos_factor])

            # Accounts for wrapping effects at ra=0, calculates the distance
            # of each event to the source.
            ra_dist = np.fabs(
                (data["ra"] - source['ra'] + np.pi) % (2. * np.pi) - np.pi)
            ra_mask = ra_dist < dPhi / 2.

            spatial_mask = dec_mask & ra_mask

            veto = veto & ~spatial_mask

        return ~veto

    @staticmethod
    def assume_background(n_s, n_coincident, n_all):
        """To save time with likelihood calculation, it can be assumed that
        all events defined as "non-coincident", because of distance in space
        and time to the source, are in fact background events. This is
        equivalent to setting S=0 for all non-coincident events. IN this
        case, the likelihood can be calculated as the product of the number
        of non-coincident events, and the likelihood of an event which has S=0.

        :param n_s: Array of expected number of events
        :param n_coincident: Number of events that were not assumed to have S=0
        :param n_all: The total number of events
        :return: Log Likelihood value for the given
        """
        return (n_all - n_coincident) * np.log1p(-n_s / n_all)

    def create_kwargs(self, data):
        kwargs = dict()
        return kwargs

    def create_llh_function(self, data):
        """Creates a likelihood function to minimise, based on the dataset.

        :param data: Dataset
        :return: LLH function that can be minimised
        """

        kwargs = self.create_kwargs(data)

        def test_statistic(params, weights):
            return self.calculate_test_statistic(
                params, weights, **kwargs)

        return test_statistic

    def calculate_test_statistic(self, params, weights, **kwargs):
        pass

    @staticmethod
    def return_llh_parameters(llh_dict):
        seeds = []
        bounds = []
        names = []

        return seeds, names, bounds

    @staticmethod
    def return_injected_parameters(mh_dict):
        return {}


@LLH.register_subclass('spatial')
class SpatialLLH(LLH):
    """Most basic LLH, in which only spatial, and otionally slso temporal,
    information is included. No Energy PDF is used, and no energy weighting
    is applied.

    """
    fit_energy = False

    def __init__(self, season, sources, llh_dict):
        LLH.__init__(self, season, sources, llh_dict)

        if "LLH Energy PDF" in llh_dict.keys():
            raise Exception("Found 'LLH Energy PDF' entry in llh_dict, "
                            "but SpatialLLH does not use Energy PDFs. \n"
                            "Please remove this entry, and try again.")

    def create_acceptance_function(self):
        """In the most simple case of spatial-only weighting, you would
        neglect the energy weighting of events. Then, you can simply assume
        that the detector acceptance is roughly proportional to the data rate,
        i.e assuming that the incident background atmospheric neutrino flux
        is uniform. Thus the acceptance of the detector is simply the
        background spatial PDF (which is a spline fitted to data as a
        function of declination). This method does, admittedly neglect the
        fact that background in the southern hemisphere is mainly composed of
        muon bundles, rather than atmospheric neutrinos. Still, it's slighty
        better than assuming a uniform detector acceptance

        :return: 1D linear interpolation
        """
        exp = data_loader(self.season["exp_path"])
        data_rate = float(len(exp))
        del exp

        # return lambda x: data_rate
        return lambda x: np.exp(self.bkg_spatial(np.sin(x))) * data_rate

    def create_llh_function(self, data):
        """Creates a likelihood function to minimise, based on the dataset.

        :param data: Dataset
        :return: LLH function that can be minimised
        """
        n_all = float(len(data))
        SoB_spacetime = []

        assumed_bkg_mask = np.ones(len(data), dtype=np.bool)

        for i, source in enumerate(np.sort(self.sources,
                                           order="Distance (Mpc)")):

            s_mask = self.select_spatially_coincident_data(data, [source])

            assumed_bkg_mask *= ~s_mask
            coincident_data = data[s_mask]

            if len(coincident_data) > 0:

                sig = self.signal_pdf(source, coincident_data)

                bkg = np.array(self.background_pdf(source, coincident_data))

                SoB_spacetime.append(sig/bkg)
                del sig
                del bkg

        n_coincident = np.sum(~assumed_bkg_mask)

        SoB_spacetime = np.array(SoB_spacetime)

        def test_statistic(params, weights):

            return self.calculate_test_statistic(
                params, weights, n_all=n_all, n_coincident=n_coincident,
                SoB_spacetime=SoB_spacetime)

        return test_statistic

    def calculate_test_statistic(self, params, weights, **kwargs):
        """Calculates the test statistic, given the parameters. Uses numexpr
        for faster calculations.

        :param params: Parameters from Minimisation
        :param weights: Normalised fraction of n_s allocated to each source
        :return: 2 * llh value (Equal to Test Statistic)
        """

        n_s = np.array(params)

        # Calculates the expected number of signal events for each source in
        # the season
        all_n_j = (n_s * weights.T[0])

        x = []

        for i, n_j in enumerate(all_n_j):
            x.append(1 + ((n_j / kwargs["n_all"]) *
                          (kwargs["SoB_spacetime"][i] - 1.)))

        if np.sum([np.sum(x_row <= 0.) for x_row in x]) > 0:
            llh_value = -50. + all_n_j

        else:

            llh_value = np.array([np.sum(np.log(y)) for y in x])

            llh_value += self.assume_background(
                all_n_j, kwargs["n_coincident"], kwargs["n_all"])

            if np.logical_and(np.sum(all_n_j) < 0,
                              np.sum(llh_value) < np.sum(-50. + all_n_j)):
                llh_value = -50. + all_n_j

        # Definition of test statistic
        return 2. * np.sum(llh_value)


@LLH.register_subclass('fixed_energy')
class FixedEnergyLLH(LLH):

    fit_energy = False

    def __init__(self, season, sources, llh_dict):

        try:
            e_pdf_dict = llh_dict["LLH Energy PDF"]
            self.energy_pdf = EnergyPDF.create(e_pdf_dict)
        except KeyError:
            raise KeyError("LLH with energy term selected, but no energy PDF "
                           "has been provided. Please add an 'LLH Energy "
                           "PDF' dictionary to the LLH dictionary, and try "
                           "again.")

        LLH.__init__(self, season, sources, llh_dict)

    def create_acceptance_function(self):
        print "Building acceptance functions in sin(dec) bins " \
              "(with fixed energy weighting)"

        mc = data_loader(self.season["mc_path"])

        acc = np.ones_like(dec_range, dtype=np.float)

        for i, dec in enumerate(dec_range):

            # Sets half width of band
            dec_width = np.deg2rad(5.)

            # Sets a declination band 5 degrees above and below the source
            min_dec = max(-np.pi / 2., dec - dec_width)
            max_dec = min(np.pi / 2., dec + dec_width)
            # Gives the solid angle coverage of the sky for the band
            omega = 2. * np.pi * (np.sin(max_dec) - np.sin(min_dec))

            band_mask = np.logical_and(np.greater(mc["trueDec"], min_dec),
                                       np.less(mc["trueDec"], max_dec))

            cut_mc = mc[band_mask]
            weights = self.energy_pdf.weight_mc(cut_mc)
            acc[i] = np.sum(weights / omega)

        f = scipy.interpolate.interp1d(
            dec_range, acc, kind='linear')

        del mc

        return f

    def create_kwargs(self, data):
        """Creates a likelihood function to minimise, based on the dataset.

        :param data: Dataset
        :return: LLH function that can be minimised
        """
        kwargs = dict()
        kwargs["n_all"] = float(len(data))
        SoB = []

        assumed_bkg_mask = np.ones(len(data), dtype=np.bool)

        ratio_spline = create_2d_ratio_spline(
            exp=data_loader(self.season["exp_path"]),
            mc=data_loader(self.season["mc_path"]),
            sin_dec_bins=self.season["sinDec bins"],
            weight_function=self.energy_pdf.weight_mc
        )

        for i, source in enumerate(np.sort(self.sources,
                                           order="Distance (Mpc)")):

            s_mask = self.select_spatially_coincident_data(data, [source])

            assumed_bkg_mask *= ~s_mask
            coincident_data = data[s_mask]

            if len(coincident_data) > 0:

                sig = self.signal_pdf(source, coincident_data)
                bkg = np.array(self.background_pdf(source, coincident_data))

                SoB_energy_ratio = [
                    np.exp(ratio_spline(x["logE"], x["sinDec"]))[0][0]
                    for x in coincident_data
                ]

                SoB.append(SoB_energy_ratio * sig/bkg)

        kwargs["n_coincident"] = np.sum(~assumed_bkg_mask)

        kwargs["SoB"] = np.array(SoB)
        return kwargs


        # def test_statistic(params, weights):
        #
        #     return self.calculate_test_statistic(
        #         params, weights, n_all=n_all, n_coincident=n_coincident,
        #         SoB=SoB)
        #
        # return test_statistic

    def calculate_test_statistic(self, params, weights, **kwargs):
        """Calculates the test statistic, given the parameters. Uses numexpr
        for faster calculations.

        :param params: Parameters from Minimisation
        :param weights: Normalised fraction of n_s allocated to each source
        :return: 2 * llh value (Equal to Test Statistic)
        """
        n_s = np.array(params)

        # Calculates the expected number of signal events for each source in
        # the season
        all_n_j = (n_s * weights.T[0])

        x = []

        for i, n_j in enumerate(all_n_j):
            x.append(1 + ((n_j / kwargs["n_all"]) *
                          (kwargs["SoB"][i] - 1.)))

        if np.sum([np.sum(x_row <= 0.) for x_row in x]) > 0:
            llh_value = -50. + all_n_j

        else:

            llh_value = np.array([np.sum(np.log(y)) for y in x])

            llh_value += self.assume_background(
                all_n_j, kwargs["n_coincident"], kwargs["n_all"])

            if np.logical_and(np.sum(all_n_j) < 0,
                              np.sum(llh_value) < np.sum(-50. + all_n_j)):
                llh_value = -50. + all_n_j

        # Definition of test statistic
        return 2. * np.sum(llh_value)


@LLH.register_subclass('standard')
class StandardLLH(FixedEnergyLLH):

    fit_energy = True

    def __init__(self, season, sources, llh_dict):
        FixedEnergyLLH.__init__(self, season, sources, llh_dict)

        # Bins for energy Log(E/GeV)
        self.energy_bins = np.linspace(1., 10., 40 + 1)

        # Sets precision for energy SoB
        self.precision = .1

        print "Loading Log(Signal/Background) Splines."

        self.SoB_spline_2Ds = load_spline(self.season)

        # print "Loaded", len(self.SoB_spline_2Ds), "Splines."

        self.acceptance_f = self.create_acceptance_function()

    def _around(self, value):
        """Produces an array in which the precision of the value
        is rounded to the nearest integer. This is then multiplied
        by the precision, and the new value is returned.

        :param value: value to be processed
        :return: value after processed
        """
        return np.around(float(value) / self.precision) * self.precision

    def create_acceptance_function(self):
        """Creates a 2D linear interpolation of the acceptance of the detector
        for the given season, as a function of declination and gamma. Returns
        this interpolation function.

        :return: 2D linear interpolation
        """

        acc_path = acceptance_path(self.season)

        with open(acc_path) as f:
            acc_dict = Pickle.load(f)

        dec_bins = acc_dict["dec"]
        gamma_bins = acc_dict["gamma"]
        values = acc_dict["acceptance"]
        f = scipy.interpolate.interp2d(
            dec_bins, gamma_bins, values.T, kind='linear')
        return f

    def acceptance(self, source, params=None):
        """Calculates the detector acceptance for a given source, using the
        2D interpolation of the acceptance as a function of declination and
        gamma. If gamma IS NOT being fit, uses the default value of gamma for
        weighting (determined in __init__). If gamma IS being fit, it will be
        the last entry in the parameter array, and is the acceptance uses
        this value.

        :param source: Source to be considered
        :param params: Parameter array
        :return: Value for the acceptance of the detector, in the given
        season, for the source
        """
        dec = source["dec"]
        gamma = params[-1]

        return self.acceptance_f(dec, gamma)

    def create_kwargs(self, data):

        kwargs = dict()

        kwargs["n_all"] = float(len(data))
        SoB_spacetime = []
        SoB_energy_cache = []

        assumed_background_mask = np.ones(len(data), dtype=np.bool)

        for i, source in enumerate(np.sort(self.sources,
                                           order="Distance (Mpc)")):

            s_mask = self.select_spatially_coincident_data(data, [source])

            assumed_background_mask *= ~s_mask
            coincident_data = data[s_mask]

            if len(coincident_data) > 0:

                sig = self.signal_pdf(source, coincident_data)
                bkg = np.array(self.background_pdf(source, coincident_data))

                SoB_spacetime.append(sig/bkg)
                del sig
                del bkg

                energy_cache = self.create_SoB_energy_cache(coincident_data)

                SoB_energy_cache.append(energy_cache)

            # print n_bkg

        kwargs["n_coincident"] = np.sum(~assumed_background_mask)

        kwargs["SoB_spacetime"] = np.array(SoB_spacetime)
        kwargs["SoB_energy_cache"] = SoB_energy_cache

        return kwargs

    def calculate_test_statistic(self, params, weights, **kwargs):
        """Calculates the test statistic, given the parameters. Uses numexpr
        for faster calculations.

        :param params: Parameters from Minimisation
        :param weights: Normalised fraction of n_s allocated to each source
        :return: 2 * llh value (Equal to Test Statistic)
        """
        n_s = np.array(params[:-1])
        gamma = params[-1]
        SoB_energy = np.array([self.estimate_energy_weights(gamma, x)
                               for x in kwargs["SoB_energy_cache"]])

        # Calculates the expected number of signal events for each source in
        # the season
        all_n_j = (n_s * weights.T[0])

        x = []

        # If n_s if negative, then removes the energy term from the likelihood

        for i, n_j in enumerate(all_n_j):
            # Switches off Energy term for negative n_s, which should in theory
            # be a continuous change that does not alter the likelihood for
            # n_s > 0 (as it is not included for n_s=0).
            if n_j < 0:
                x.append(1 + ((n_j / kwargs["n_all"]) * (
                        kwargs["SoB_spacetime"][i] - 1.)))
            else:

                x.append(1 + ((n_j / kwargs["n_all"]) * (
                    SoB_energy[i] * kwargs["SoB_spacetime"][i] - 1.)))


        if np.sum([np.sum(x_row <= 0.) for x_row in x]) > 0:
            llh_value = -50. + all_n_j

        else:

            llh_value = np.array([np.sum(np.log(y)) for y in x])

            # print "llh value", llh_value, n_s, kwargs["n_all"]
            # raw_input("prompt")

            llh_value += self.assume_background(
                all_n_j, kwargs["n_coincident"], kwargs["n_all"])

            if np.logical_and(np.sum(all_n_j) < 0,
                              np.sum(llh_value) < np.sum(-50. + all_n_j)):
                llh_value = -50. + all_n_j

        # Definition of test statistic
        return 2. * np.sum(llh_value)


# ==============================================================================
# Energy Log(Signal/Background) Ratio
# ==============================================================================

    def create_SoB_energy_cache(self, cut_data):
        """Evaluates the Log(Signal/Background) values for all coincident
        data. For each value of gamma in self.gamma_support_points, calculates
        the Log(Signal/Background) values for the coincident data. Then saves
        each weight array to a dictionary.

        :param cut_data: Subset of the data containing only coincident events
        :return: Dictionary containing SoB values for each event for each
        gamma value.
        """
        energy_SoB_cache = dict()

        for gamma in self.SoB_spline_2Ds.keys():
            energy_SoB_cache[gamma] = self.SoB_spline_2Ds[gamma].ev(
                cut_data["logE"], cut_data["sinDec"])

        return energy_SoB_cache

    def estimate_energy_weights(self, gamma, energy_SoB_cache):
        """Quickly estimates the value of Signal/Background for Gamma.
        Uses pre-calculated values for first and second derivatives.
        Uses a Taylor series to estimate S(gamma), unless SoB has already
        been calculated for a given gamma.

        :param gamma: Spectral Index
        :param energy_SoB_cache: Weight cache
        :return: Estimated value for S(gamma)
        """
        if gamma in energy_SoB_cache.keys():
            val = np.exp(energy_SoB_cache[gamma])
        else:
            g1 = self._around(gamma)
            dg = self.precision

            g0 = self._around(g1 - dg)
            g2 = self._around(g1 + dg)

            # Uses Numexpr to quickly estimate S(gamma)

            S0 = energy_SoB_cache[g0]
            S1 = energy_SoB_cache[g1]
            S2 = energy_SoB_cache[g2]

            val = numexpr.evaluate(
                "exp((S0 - 2.*S1 + S2) / (2. * dg**2) * (gamma - g1)**2" + \
                " + (S2 -S0) / (2. * dg) * (gamma - g1) + S1)"
            )

        return val

    @staticmethod
    def return_llh_parameters(llh_dict):
        e_pdf = EnergyPDF.create(llh_dict["LLH Energy PDF"])
        return e_pdf.return_energy_parameters()

    @staticmethod
    def return_injected_parameters(mh_dict):

        try:

            inj = mh_dict["inj kwargs"]["Injection Energy PDF"]
            llh = mh_dict["llh_dict"]["LLH Energy PDF"]

            if inj["Name"] == llh["Name"]:
                e_pdf = EnergyPDF.create(inj)
                return e_pdf.return_injected_parameters()

        except KeyError:
            pass

        seeds, bounds, names = LLH.get_parameters(mh_dict["llh_dict"])

        res_dict = {}
        for key in names:
            res_dict[key] = np.nan

        return res_dict


def generate_dynamic_flare_class(season, sources, llh_dict):

    try:
        mh_name = llh_dict["name"]
    except KeyError:
        raise KeyError("No LLH specified.")

    # Set up dynamic inheritance

    try:
        ParentLLH = LLH.subclasses[mh_name]
    except KeyError:
        raise KeyError("Parent class {} not found.".format(mh_name))

    # Defines custom Flare class

    class FlareLLH(ParentLLH):

        def create_flare_llh_function(self, data, flare_veto,
                                      n_all, src, n_season):

            coincident_data = data[~flare_veto]
            kwargs = self.create_kwargs(coincident_data)
            kwargs["n_all"] = n_all
            weights = np.array([1.])

            def test_statistic(params):
                return self.calculate_test_statistic(
                    params, weights, **kwargs)



            # def base_f(params):
            #     return test_statistic(coincident_data)(
            #         params)

            # Super ugly-looking code that magically takes the old llh
            # object, sets the assume_background contribution to zero,
            # and then adds on a new assume_season_background where mutiple
            # datasets are treated as one season of data.

            def combined_test_statistic(params):
                return test_statistic(params) + (
                        2 * self.assume_season_background(
                    params[0], np.sum(~flare_veto), n_season, n_all)
                )

            # base_ts =
            #
            # base_ts += self.assume_season_background()

            # sig = self.signal_spatial(source, coincident_data)
            # bkg = self.background_spatial(coincident_data)
            # SoB_spacetime = sig/bkg
            # del sig
            # del bkg
            #
            # # If an llh energy PDF has been provided, calculate the SoB values
            # # for the coincident data, and stores it in a cache.
            # if hasattr(self, "energy_pdf"):
            #     SoB_energy_cache = self.create_SoB_energy_cache(coincident_data)
            #
            #     # If gamma is not going to be fit, replaces the SoB energy
            #     # cache with the weight array corresponding to the gamma provided
            #     # in the llh energy PDF
            #     if not self.fit_gamma:
            #         SoB_energy_cache = self.estimate_energy_weights(
            #             self.default_gamma, SoB_energy_cache)
            #
            # else:
            #     SoB_energy_cache = None
            #
            # def test_statistic(params):
            #     return self.calculate_test_statistic(
            #         params, n_season, n_all, SoB_spacetime,
            #         SoB_energy_cache)

            return combined_test_statistic

        @staticmethod
        def assume_background(n_s, n_coincident, n_all):
            """In the standard create_llh_function method that the FlareClass
            inherits, this method will be called. To maintain modularity, we
            simply set it to 0 here. The Flare class treats all neutrino events
            collectively, rather than splitting them by season. As a result,
            the assume_season_background method is called seperately to handle
            the non-coincident neutrinos.

            :param n_s: Array of expected number of events
            :param n_coincident: Number of events that were not assumed to have S=0
            :param n_all: The total number of events
            :return: 0.
            """
            return 0.

        def signal_pdf(self, source, cut_data):
            """Calculates the value of the signal spatial PDF for a given source
            for each event in the coincident data subsample. If there is a Time PDF
            given, also calculates the value of the signal Time PDF for each event.
            Returns either the signal spatial PDF values, or the product of the
            signal spatial and time PDFs.

            :param source: Source to be considered
            :param cut_data: Subset of Dataset with coincident events
            :return: Array of Signal Spacetime PDF values
            """
            space_term = self.signal_spatial(source, cut_data)

            return space_term

        # ==========================================================================
        # Background PDF
        # ==========================================================================

        def background_pdf(self, source, cut_data):
            """For the flare search, generating repeated box time PDFs would
            be required to recalculate the
            """
            space_term = self.background_spatial(cut_data)

            return space_term

        # def calculate_test_statistic(self, params, n_all, n_coincident,
        #                              SoB_spacetime, SoB_energy_cache=None,
        #                              weights=1.):
        #     """Calculates the test statistic, given the parameters. Uses numexpr
        #     for faster calculations.
        #
        #     :param params: Parameters from minimisation
        #     :return: Test Statistic
        #     """
        #     n_mask = len(SoB_spacetime)
        #
        #     # If fitting gamma and calculates the energy weights for the given
        #     # value of gamma
        #     if self.fit_gamma:
        #         n_s = np.array(params[:-1])
        #         gamma = params[-1]
        #
        #         SoB_energy = self.estimate_energy_weights(gamma, SoB_energy_cache)
        #
        #     # If using energy information but with a fixed value of gamma,
        #     # sets the weights as equal to those for the provided gamma value.
        #     elif SoB_energy_cache is not None:
        #         n_s = np.array(params)
        #         SoB_energy = SoB_energy_cache
        #
        #     # If not using energy information, assigns a weight of 1. to each event
        #     else:
        #         n_s = np.array(params)
        #         SoB_energy = 1.
        #
        #     if len(SoB_spacetime) > 0:
        #         # Evaluate the likelihood function for neutrinos close to each source
        #         llh_value = np.sum(np.log((
        #             1 + ((n_s/n_all) * (SoB_energy * SoB_spacetime)))))
        #
        #     else:
        #         llh_value = 0.
        #
        #     llh_value += self.assume_season_background(n_s, n_mask, n_coincident, n_all)
        #
        #     # Definition of test statistic
        #     return 2. * llh_value

        @staticmethod
        def assume_season_background(n_s, n_mask, n_season, n_all):
            """To save time with likelihood calculation, it can be assumed that
            all events defined as "non-coincident", because of distance in space
            and time to the source, are in fact background events. This is
            equivalent to setting S=0 for all non-coincident events. IN this
            case, the likelihood can be calculated as the product of the number
            of non-coincident events, and the likelihood of an event which has S=0.

            :param n_s: Array of expected number of events
            :param n_mask: Number of events that were not assumed to have S=0
            :param n_all: The total number of events
            :return: Log Likelihood value for the given
            """
            return (n_season - n_mask) * np.log1p(-n_s / n_all)

        def estimate_significance(self, coincident_data, source):
            """Finds events in the coincident dataset (spatially and temporally
            overlapping sources), which are significant. This is defined as having a
            Signal/Background Ratio that is greater than 1. The S/B ratio is
            calculating using spatial and energy PDFs.

            :param coincident_data: Data overlapping the source spatially/temporally
            :param source: Source to be considered
            :return: SoB of events in coincident dataset
            """
            sig = self.signal_spatial(source, coincident_data)
            bkg = self.background_spatial(coincident_data)
            SoB_space = sig / bkg

            SoB_energy_cache = self.create_SoB_energy_cache(coincident_data)

            # ChangeME?

            SoB_energy = self.estimate_energy_weights(
                    gamma=3.0, energy_SoB_cache=SoB_energy_cache)

            SoB = SoB_space * SoB_energy
            return SoB

        def find_significant_events(self, coincident_data, source):
            """Finds events in the coincident dataset (spatially and temporally
            overlapping sources), which are significant. This is defined as having a
            Signal/Background Ratio that is greater than 1. The S/B ratio is
            calculating using spatial and energy PDFs.

            :param coincident_data: Data overlapping the source spatially/temporally
            :param source: Source to be considered
            :return: Significant events in coincident dataset
            """

            SoB = self.estimate_significance(coincident_data, source)

            mask = SoB > 1.0

            return coincident_data[mask]

    return FlareLLH(season, sources, llh_dict)


if __name__ == "__main__":
    from flarestack.shared import fs_scratch_dir
    from scipy.interpolate import InterpolatedUnivariateSpline

    g = EnergyPDF.create(
        {
            "Name": "Power Law",
            "Gamma": 2.2
        }
    )

    e_range = np.logspace(0, 7, 1e3)

    f = InterpolatedUnivariateSpline(e_range, np.log(g.f(e_range)))

    path = fs_scratch_dir + "tester_spline.npy"

    print path

    with open(path, "wb") as h:
        Pickle.dump(f, h)

    e_pdf = {
        "Name": "Spline",
        "Spline Path": path,
    }

    from flarestack.data.icecube.ps_tracks.ps_v002_p01 import IC86_1_dict
    from flarestack.utils.prepare_catalogue import ps_catalogue_name
    from flarestack.core.injector import MockUnblindedInjector, Injector

    llh_dict = {
        "Name": "FixedEnergy",
        "LLH Time PDF": {
            "Name": "Steady"
        },
        "LLH Energy PDF": {
            "Name": "Power Law",
            "Gamma": 2.0
        }
        # "LLH Energy PDF": e_pdf_dict
    }
    source = np.load(ps_catalogue_name(0.0))

    llh = LLH.create(IC86_1_dict, source, llh_dict)

    # inj = MockUnblindedInjector(IC86_1_dict, source)
    inj = Injector(IC86_1_dict, source)

    data = inj.create_dataset(0.0)
    f = llh.create_llh_function(data)

    weights = np.array([1.])

    for i in np.linspace(0.0, 10.0, 21):
        print i, f([i], weights)

