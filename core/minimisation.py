import numpy as np
import resource
import random
import os, os.path
import argparse
import cPickle as Pickle
import scipy.optimize
from core.injector import Injector, MockUnblindedInjector
from core.llh import LLH, FlareLLH
from shared import name_pickle_output_dir, fit_setup, inj_dir_name,\
    plot_output_dir, scale_shortener
import matplotlib.pyplot as plt


def time_smear(inj):
    inj_time = inj["Injection Time PDF"]
    max_length = inj_time["Max Offset"] - inj_time["Min Offset"]
    offset = np.random.random() * max_length + inj_time["Min Offset"]
    inj_time["Offset"] = offset
    return inj_time


class MinimisationHandler:
    """Generic Class to handle both dataset creation and llh minimisation from
    experimental data and Monte Carlo simulation. Initialised with a set of
    IceCube datasets, a list of sources, and independent sets of arguments for
    the injector and the likelihood.
    """
    n_trials_default = 1000

    def __init__(self, mh_dict):

        sources = np.sort(np.load(mh_dict["catalogue"]), order="Distance (Mpc)")

        self.name = mh_dict["name"]
        self.pickle_output_dir = name_pickle_output_dir(self.name)
        self.injectors = dict()
        self.llhs = dict()
        self.seasons = mh_dict["datasets"]
        self.sources = sources

        # Checks whether signal injection should be done with a sliding PDF
        # within a larger window, or remain fixed at the specified time

        inj = dict(mh_dict["inj kwargs"])

        try:
            self.time_smear = inj["Injection Time PDF"]["Time Smear?"]
        except KeyError:
            self.time_smear = False

        if self.time_smear:
            inj["Injection Time PDF"] = time_smear(inj)

        self.inj_kwargs = inj
        self.llh_kwargs = mh_dict["llh kwargs"]

        # Checks if the code should search for flares. By default, this is
        # not done.
        try:
            self.flare = self.llh_kwargs["Flare Search?"]
        except KeyError:
            self.flare = False

        if self.flare:
            self.run = self.run_flare
        else:
            self.run = self.run_stacked

        # Checks if minimiser should be seeded from a brute scan

        try:
            self.brute = self.llh_kwargs["Brute Seed?"]
        except KeyError:
            self.brute = False

        # Checks if source weights should be fitted individually

        try:
            self.fit_weights = self.llh_kwargs["Fit Weights?"]

            if "Negative n_s?" in self.llh_kwargs.keys():
                if self.llh_kwargs["Negative n_s?"]:
                    raise Exception("Attempted to mix fitting weights with "
                                    "negative n_s. The code is not able to "
                                    "handle this setup!")

        except KeyError:
            self.fit_weights = False

        if self.fit_weights:
            self.run_trial = self.run_fit_weight_trial
        else:
            self.run_trial = self.run_fixed_weight_trial

        # Checks if data should be "mock-unblinded", where a fixed seed
        # background scramble is done for injection stage. Only calls to the
        # Unblinder class will have this attribute.

        if hasattr(self, "unblind_dict"):
            self.mock_unblind = mh_dict["Mock Unblind"]
        else:
            self.mock_unblind = False

        # For each season, we create an independent injector and a
        # likelihood, using the source list along with the sets of energy/time
        # PDFs provided in inj_kwargs and llh_kwargs.
        for season in self.seasons:

            if not self.flare:
                self.llhs[season["Name"]] = LLH(season, sources,
                                                **self.llh_kwargs)
            else:
                self.llhs[season["Name"]] = FlareLLH(season, sources,
                                                     **self.llh_kwargs)

            if self.mock_unblind:
                self.injectors[season["Name"]] = MockUnblindedInjector(
                    season, sources, **self.inj_kwargs)
            else:
                self.injectors[season["Name"]] = Injector(
                    season, sources, **self.inj_kwargs)

        p0, bounds, names = fit_setup(self.llh_kwargs, sources, self.flare)

        self.p0 = p0
        self.bounds = bounds
        self.param_names = names

        # Sets the default flux scale for finding sensitivity
        # Default value is 1 (Gev)^-1 (cm)^-2 (s)^-1
        self.bkg_ts = None

        # self.clean_true_param_values()

    def clear(self):

        self.injectors.clear()
        self.llhs.clear()

        del self

    def dump_results(self, results, scale, seed):
        """Takes the results of a set of trials, and saves the dictionary as
        a pickle pkl_file. The flux scale is used as a parent directory, and the
        pickle pkl_file itself is saved with a name equal to its random seed.

        :param results: Dictionary of Minimisation results from trials
        :param scale: Scale of inputted flux
        :param seed: Random seed used for running of trials
        """

        write_dir = self.pickle_output_dir + scale_shortener(scale) + "/"

        # Tries to create the parent directory, unless it already exists
        try:
            os.makedirs(write_dir)
        except OSError:
            pass

        file_name = write_dir + str(seed) + ".pkl"

        print "Saving to", file_name

        with open(file_name, "wb") as f:
            Pickle.dump(results, f)

    def dump_injection_values(self, scale):

        inj_dict = dict()
        for source in self.sources:
            name = source["Name"]
            n_inj = 0
            for inj in self.injectors.itervalues():
                try:
                    n_inj += inj.ref_fluxes[scale_shortener(scale)][name]

                # If source not overlapping season, will not be in dict
                except KeyError:
                    pass

            default = {
                "n_s": n_inj
            }

            if "Gamma" in self.param_names:
                default["Gamma"] = self.inj_kwargs["Injection Energy PDF"][
                    "Gamma"]

            if self.flare:
                fs = [inj.time_pdf.sig_t0(source)
                      for inj in self.injectors.itervalues()]
                true_fs = min(fs)
                fe = [inj.time_pdf.sig_t1(source)
                      for inj in self.injectors.itervalues()]
                true_fe = max(fe)

                if self.time_smear:
                    inj_time = self.inj_kwargs["Injection Time PDF"]
                    offset = inj_time["Offset"]
                    true_fs -= offset
                    true_fe -= offset

                    min_offset = inj_time["Min Offset"]
                    max_offset = inj_time["Max Offset"]
                    med_offset = 0.5*(max_offset + min_offset)

                    true_fs += med_offset
                    true_fe += med_offset

                true_l = true_fe - true_fs

                sim = [
                    list(np.random.uniform(true_fs, true_fe,
                                      np.random.poisson(n_inj)))
                    for _ in range(1000)
                ]

                s = []
                e = []
                l = []

                for data in sim:
                    if data != []:
                        s.append(min(data))
                        e.append(max(data))
                        l.append(max(data) - min(data))

                if len(s) > 0:
                    med_s = np.median(s)
                    med_e = np.median(e)
                    med_l = np.median(l)
                else:
                    med_s = np.nan
                    med_e = np.nan
                    med_l = np.nan

                print med_s, med_e, med_l

                default["Flare Start"] = med_s
                default["Flare End"] = med_e
                default["Flare Length"] = med_l

            inj_dict[name] = default

        inj_dir = inj_dir_name(self.name)

        # Tries to create the parent directory, unless it already exists
        try:
            os.makedirs(inj_dir)
        except OSError:
            pass

        file_name = inj_dir + scale_shortener(scale) + ".pkl"
        with open(file_name, "wb") as f:
            Pickle.dump(inj_dict, f)

    def iterate_run(self, scale=1, n_steps=5, n_trials=50):

        scale_range = np.linspace(0., scale, n_steps)[1:]

        self.run(n_trials*10, scale=0.0)

        for scale in scale_range:
            self.run(n_trials, scale)

    def run_stacked(self, n_trials=n_trials_default, scale=1.):

        seed = int(random.random() * 10 ** 8)
        np.random.seed(seed)

        param_vals = [[] for x in self.p0]
        ts_vals = []
        flags = []

        print "Generating", n_trials, "trials!"

        for i in range(int(n_trials)):
            raw_f = self.run_trial(scale)

            def llh_f(scale):
                return -np.sum(raw_f(scale))

            if self.brute:

                brute_range = [
                    (max(x, -30), min(y, 30)) for (x, y) in self.bounds]

                start_seed = scipy.optimize.brute(
                    llh_f, ranges=brute_range, finish=None, Ns=40)
            else:
                start_seed = self.p0

            res = scipy.optimize.minimize(
                llh_f, start_seed, bounds=self.bounds)

            flag = res.status
            vals = res.x
            # ts = -res.fun
            # If the minimiser does not converge, repeat with brute force
            if flag == 1:

                vals = scipy.optimize.brute(llh_f, ranges=self.bounds,
                                            finish=None)
                # ts = -llh_f(vals)

            if res is not None:

                for j, val in enumerate(list(vals)):
                    param_vals[j].append(val)

                best_llh = raw_f(vals)

                if self.fit_weights:

                    ts = np.sum(best_llh)

                else:
                    ts = np.sign(vals[0]) * np.sum(best_llh)
                    # print best_llh, vals[0], ts
                    # raw_input("prompt")

                ts_vals.append(ts)
                flags.append(flag)

        mem_use = str(
            float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1.e6)
        print ""
        print 'Memory usage max: %s (Gb)' % mem_use

        n_inj = 0
        for inj in self.injectors.itervalues():
            for val in inj.ref_fluxes[scale_shortener(scale)].itervalues():
                n_inj += val
        print ""
        print "Injected with an expectation of", n_inj, "events."

        print ""
        print "FIT RESULTS:"
        print ""

        for i, param in enumerate(param_vals):
            print "Parameter", self.param_names[i], ":", np.mean(param), \
                np.median(param), np.std(param)
        print "Test Statistic:", np.mean(ts_vals), np.median(ts_vals), np.std(
            ts_vals)
        print ""

        print "FLAG STATISTICS:"
        for i in sorted(np.unique(flags)):
            print "Flag", i, ":", flags.count(i)

        results = {
            "TS": ts_vals,
            "Parameters": param_vals,
            "Flags": flags
        }

        self.dump_results(results, scale, seed)

        self.dump_injection_values(scale)

    def make_fixed_weight_matrix(self, params):

        # Creates a matrix fixing the fraction of the total signal that
        # is expected in each Source+Season pair. The matrix is
        # normalised to 1, so that for a given total n_s, the expectation
        # for the ith season for the jth source is given by:
        #  n_exp = n_s * weight_matrix[i][j]

        src = np.sort(self.sources, order="Distance (Mpc)")
        dist_weight = src["Distance (Mpc)"] ** -2

        weights_matrix = np.ones([len(self.seasons), len(self.sources)])

        for i, season in enumerate(self.seasons):
            llh = self.llhs[season["Name"]]
            acc = []

            time_weights = []

            for source in src:
                time_weights.append(llh.time_pdf.effective_injection_time(
                    source))
                acc.append(llh.acceptance(source, params))

            acc = np.array(acc).T

            w = acc * dist_weight * np.array(time_weights)

            w = w[:, np.newaxis]

            for j, ind_w in enumerate(w.T):
                weights_matrix[i][j] = ind_w

        weights_matrix /= np.sum(weights_matrix)
        return weights_matrix

    def run_fixed_weight_trial(self, scale=1.):

        llh_functions = dict()
        n_all = dict()

        for season in self.seasons:
            dataset = self.injectors[season["Name"]].create_dataset(scale)
            llh_f = self.llhs[season["Name"]].create_llh_function(dataset)
            llh_functions[season["Name"]] = llh_f
            n_all[season["Name"]] = len(dataset)

        def f_final(raw_params):

            # If n_s is less than or equal to 0, set gamma to be 3.7 (equal to
            # atmospheric background). This is continuous at n_s=0, but fixes
            # relative weights of sources/seasons for negative n_s values.

            params = list(raw_params)

            # if (len(params) > 1) and (params[0] < 0):
            #     params[1] = 3.7

            # Calculate relative contribution of each source/season

            weights_matrix = self.make_fixed_weight_matrix(params)

            # Having created the weight matrix, loops over each season of
            # data and evaluates the TS function for that season

            ts_val = 0
            for i, season in enumerate(self.seasons):
                w = weights_matrix[i][:, np.newaxis]
                ts_val += np.sum(llh_functions[season["Name"]](params, w))

            return ts_val

        return f_final

    def run_fit_weight_trial(self, scale):
        llh_functions = dict()
        n_all = dict()

        src = np.sort(self.sources, order="Distance (Mpc)")
        dist_weight = src["Distance (Mpc)"] ** -2

        for season in self.seasons:
            dataset = self.injectors[season["Name"]].create_dataset(scale)
            llh_f = self.llhs[season["Name"]].create_llh_function(dataset)
            llh_functions[season["Name"]] = llh_f
            n_all[season["Name"]] = len(dataset)

        def f_final(params):

            # Creates a matrix fixing the fraction of the total signal that
            # is expected in each Source+Season pair. The matrix is
            # normalised to 1, so that for a given total n_s, the expectation
            # for the ith season for the jth source is given by:
            #  n_exp = n_s * weight_matrix[i][j]

            weights_matrix = np.ones([len(self.seasons), len(self.sources)])

            for i, season in enumerate(self.seasons):
                llh = self.llhs[season["Name"]]
                acc = []

                time_weights = []

                for source in src:
                    time_weights.append(
                        llh.time_pdf.effective_injection_time(
                            source))
                    acc.append(llh.acceptance(source, params))

                acc = np.array(acc).T

                w = acc * np.array(time_weights)

                w = w[:, np.newaxis]

                for j, ind_w in enumerate(w.T):
                    weights_matrix[i][j] = ind_w

            for i, row in enumerate(weights_matrix.T):
                if np.sum(row) > 0:
                    row /= np.sum(row)

            # weights_matrix /= np.sum(weights_matrix)

            # Having created the weight matrix, loops over each season of
            # data and evaluates the TS function for that season

            ts_val = 0
            for i, season in enumerate(self.seasons):
                w = weights_matrix[i][:, np.newaxis]
                ts_val += llh_functions[season["Name"]](params, w)

            return ts_val

        return f_final

    def run_flare(self, n_trials=n_trials_default, scale=1.):
        """Runs iterations of a flare search, and dumps results as pickle files.
        For stacking of multiple soyrces, due to computational constraints,
        each flare is minimised entirely independently. The TS values from
        each flare is then summed to give an overall TS value. The results
        for each source are stored separately.

        :param n_trials: Number of trials to perform
        :param scale: Flux scale
        """

        # Selects the key corresponding to time for the given IceCube dataset
        # (enables use of different data samples)

        time_key = self.seasons[0]["MJD Time Key"]

        seed = int(random.random() * 10 ** 8)
        np.random.seed(seed)

        print "Running", n_trials, "trials"

        # Initialises lists for all values that will need to be stored,
        # in order to verify that the minimisation is working successfuly

        results = {
            "TS": []
        }

        for source in self.sources:
            results[source["Name"]] = {
                "TS": [],
                "Parameters": [[] for x in self.param_names]
            }

        # Loop over trials

        for _ in range(int(n_trials)):

            datasets = dict()

            full_data = dict()

            # Loop over each data season

            for season in self.seasons:

                # Generate a scrambled dataset, and save it to the datasets
                # dictionary. Loads the llh for the season.

                data = self.injectors[season["Name"]].create_dataset(scale)
                llh = self.llhs[season["Name"]]

                full_data[season["Name"]] = data

                # Loops over each source in catalogue

                for source in self.sources:

                    # Identify spatially- and temporally-coincident data

                    mask = llh.select_spatially_coincident_data(data, [source])
                    spatial_coincident_data = data[mask]

                    t_mask = np.logical_and(
                        np.greater(
                            spatial_coincident_data[time_key],
                            llh.time_pdf.sig_t0(source)),
                        np.less(
                            spatial_coincident_data[time_key],
                            llh.time_pdf.sig_t1(source))
                    )

                    coincident_data = spatial_coincident_data[t_mask]

                    # Creates empty dictionary to save info

                    name = source["Name"]
                    if name not in datasets.keys():
                        datasets[name] = dict()

                    # If there are events in the window...

                    if len(coincident_data) > 0:

                        new_entry = dict(season)
                        new_entry["Coincident Data"] = coincident_data
                        new_entry["Start (MJD)"] = llh.time_pdf.t0
                        new_entry["End (MJD)"] = llh.time_pdf.t1

                        # Identify significant events (S/B > 1)

                        significant = llh.find_significant_events(
                            coincident_data, source)

                        new_entry["Significant Times"] = significant[time_key]

                        new_entry["N_all"] = len(data)

                        datasets[name][season["Name"]] = new_entry

            stacked_ts = 0.0

            # Minimisation of each source

            for (source, source_dict) in datasets.iteritems():

                src = self.sources[self.sources["Name"] == source]

                # Create a full list of all significant times

                all_times = []
                n_tot = 0
                for season_dict in source_dict.itervalues():
                    new_times = season_dict["Significant Times"]
                    all_times.extend(new_times)
                    n_tot += len(season_dict["Coincident Data"])

                all_times = np.array(sorted(all_times))

                # Minimum flare duration (days)
                min_flare = 0.25
                # Conversion to seconds
                min_flare *= 60 * 60 * 24

                # Length of search window in livetime

                search_window = np.sum([
                    llh.time_pdf.effective_injection_time(src)
                    for llh in self.llhs.itervalues()]
                )

                # If a maximum flare length is specified, sets that here

                if "Max Flare" in self.llh_kwargs["LLH Time PDF"].keys():
                    # Maximum flare given in days, here converted to seconds
                    max_flare = self.llh_kwargs["LLH Time PDF"]["Max Flare"] * (
                        60 * 60 * 24
                    )

                else:
                    max_flare = search_window

                # Loop over all flares, and check which combinations have a
                # flare length between the maximum and minimum values

                pairs = []

                for x in all_times:
                    for y in all_times:
                        if y > x:
                            pairs.append((x, y))

                # If there is are no pairs meeting this criteria, skip

                if len(pairs) == 0:
                    # print "Continuing because no pairs"
                    continue

                all_res = []
                all_ts = []

                # Loop over each possible significant neutrino pair

                for pair in pairs:
                    t_start = pair[0]
                    t_end = pair[1]

                    # Calculate the length of the neutrino flare in livetime

                    flare_time = np.array(
                        (t_start, t_end),
                        dtype=[
                            ("Start Time (MJD)", np.float),
                            ("End Time (MJD)", np.float),
                        ]
                    )

                    flare_length = np.sum([
                        llh.time_pdf.effective_injection_time(flare_time)
                        for llh in self.llhs.itervalues()]
                    )

                    # If the flare is between the minimum and maximum length

                    if flare_length < min_flare:
                        # print "Continuing because flare too short"
                        continue
                    elif flare_length > max_flare:
                        # print "Continuing because flare too long"
                        continue

                    # Marginalisation term is length of flare in livetime
                    # divided by max flare length in livetime. Accounts
                    # for the additional short flares that can be fitted
                    # into a given window

                    overall_marginalisation = flare_length / max_flare

                    # Each flare is evaluated accounting for the
                    # background on the sky (the non-coincident
                    # data), which is given by the number of
                    # neutrinos on the sky during the given
                    # flare. (NOTE THAT IT IS NOT EQUAL TO THE
                    # NUMBER OF NEUTRINOS IN THE SKY OVER THE
                    # ENTIRE SEARCH WINDOW)

                    n_all = np.sum([np.sum(~np.logical_or(
                        np.less(data[time_key], t_start),
                        np.greater(data[time_key], t_end)))
                                    for data in full_data.itervalues()])

                    llhs = dict()

                    # Loop over data seasons

                    for i, (name, season_dict) in enumerate(
                            sorted(source_dict.iteritems())):

                        llh = self.llhs[season_dict["Name"]]

                        # Check that flare overlaps with season

                        inj_time = llh.time_pdf.effective_injection_time(
                            flare_time
                        )

                        if not inj_time > 0:
                            # print "Continuing due to season having 0 overlap"
                            continue

                        coincident_data = season_dict["Coincident Data"]

                        data = full_data[name]

                        n_season = np.sum(~np.logical_or(
                            np.less(data[time_key], t_start),
                            np.greater(data[time_key], t_end)))

                        # Removes non-coincident data

                        flare_veto = np.logical_or(
                            np.less(coincident_data[time_key], t_start),
                            np.greater(coincident_data[time_key], t_end)
                        )

                        # Sets start and end time of overlap of flare
                        # with given season.

                        t_s = min(
                            max(t_start, season_dict["Start (MJD)"]),
                            season_dict["End (MJD)"])
                        t_e = max(
                            min(t_end, season_dict["End (MJD)"]),
                            season_dict["Start (MJD)"])

                        # Checks to make sure that there are
                        # neutrinos in the sky at all. There should
                        # be, due to the definition of the flare window.

                        if n_all > 0:
                            pass
                        else:
                            raise Exception("Events are leaking "
                                            "somehow!")

                        # Creates the likelihood function for the flare

                        flare_f = llh.create_flare_llh_function(
                            coincident_data, flare_veto, n_all, src, n_season)

                        llhs[season_dict["Name"]] = {
                            # "llh": unblind_llh,
                            "f": flare_f,
                            "flare length": flare_length
                        }

                    # From here, we have normal minimisation behaviour

                    def f_final(params):

                        # weights_matrix = np.ones(len(llhs))

                        # for i, (name, llh_dict) in enumerate(
                        #         sorted(llhs.iteritems())):
                        #     T = llh_dict["flare length"]
                        #     acc = llh.acceptance(src, params)
                        #     weights_matrix[i] = T * acc
                        #
                        # weights_matrix /= np.sum(weights_matrix)

                        # Marginalisation is done once, not per-season

                        ts = 2 * np.log(overall_marginalisation)

                        for i, (name, llh_dict) in enumerate(
                                sorted(llhs.iteritems())):

                            ts += llh_dict["f"](params)

                        return -ts

                    res = scipy.optimize.fmin_l_bfgs_b(
                        f_final, self.p0, bounds=self.bounds,
                        approx_grad=True)

                    all_res.append(res)
                    all_ts.append(-res[1])

                max_ts = max(all_ts)
                stacked_ts += max_ts
                index = all_ts.index(max_ts)

                best_start = pairs[index][0]
                best_end = pairs[index][1]

                best_time = np.array(
                    (best_start, best_end),
                    dtype=[
                        ("Start Time (MJD)", np.float),
                        ("End Time (MJD)", np.float),
                    ]
                )

                best_length = np.sum([
                    llh.time_pdf.effective_injection_time(best_time)
                    for llh in self.llhs.itervalues()]
                ) / (60 * 60 * 24)

                best = [x for x in all_res[index][0]] + [
                    best_start, best_end, best_length
                ]

                for k, val in enumerate(best):
                    results[source]["Parameters"][k].append(val)

                results[source]["TS"].append(max_ts)

            results["TS"].append(stacked_ts)

        mem_use = str(
            float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1.e6)
        print ""
        print 'Memory usage max: %s (Gb)' % mem_use

        full_ts = results["TS"]

        print "Combined Test Statistic:"
        print np.mean(full_ts), np.median(full_ts), np.std(
              full_ts)

        for source in self.sources:
            print "Results for", source["Name"]

            combined_res = results[source["Name"]]

            full_ts = combined_res["TS"]

            full_params = np.array(combined_res["Parameters"])

            for i, column in enumerate(full_params):
                print self.param_names[i], ":", np.mean(column),\
                    np.median(column), np.std(column)

            print "Test Statistic", np.mean(full_ts), np.median(full_ts), \
                np.std(full_ts), "\n"

        self.dump_results(results, scale, seed)

        self.dump_injection_values(scale)

    def scan_likelihood(self, scale=1.):
        """Generic wrapper to perform a likelihood scan a background scramble
        with an injection of signal given by scale.

        :param scale: Flux scale to inject
        """
        f = self.run_trial(scale)

        def g(x):
            return -np.sum(f(x))
        #
        # g = self.run_trial(scale)

        res = scipy.optimize.minimize(
            g, self.p0, bounds=self.bounds)

        print "Scan results:"
        print res
        #
        # raw_input("prompt")

        plt.figure(figsize=(8, 4 + 2*len(self.p0)))

        for i, bounds in enumerate(self.bounds):
            plt.subplot(len(self.p0), 1, 1 + i)

            best = list(res.x)

            n_range = np.linspace(max(bounds[0], -100),
                                  min(bounds[1], 100), 1e2)
            y = []

            for n in n_range:

                best[i] = n

                new = g(best)
                try:
                    y.append(new[0][0])
                except IndexError:
                    y.append(new)

            plt.plot(n_range, y)
            plt.xlabel(self.param_names[i])
            # plt.ylabel(r"$-\log(\frac{\mathcal{L}}{\mathcal{L_{0}})$")
            plt.ylabel(r"$-\log(\mathcal{L}/\mathcal{L}_{0})$")


            print "PARAM:", self.param_names[i]
            min_y = np.min(y)
            print "Minimum value of", min_y,

            min_index = y.index(min_y)
            min_n = n_range[min_index]
            print "at", min_n

            l_y = np.array(y[:min_index])
            try:
                l_y = min(l_y[l_y > (min_y + 0.5)])
                l_lim = n_range[y.index(l_y)]
            except ValueError:
                l_lim = 0

            u_y = np.array(y[min_index:])
            try:
                u_y = min(u_y[u_y > (min_y + 0.5)])
                u_lim = n_range[y.index(u_y)]
            except ValueError:
                u_lim = ">" + str(max(n_range))

            print "One Sigma interval between", l_lim, "and", u_lim

        path = plot_output_dir(self.name) + "llh_scan.pdf"

        plt.suptitle(os.path.basename(self.name[:-1]))

        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass

        plt.savefig(path)
        plt.close()

        print "Saved to", path

    def check_flare_background_rate(self):

        results = [[] for x in self.seasons]
        total = [[] for x in self.seasons]

        for i in range(int(1000)):

            # Loop over each data season

            for j, season in enumerate(sorted(self.seasons)):

                # Generate a scrambled dataset, and save it to the datasets
                # dictionary. Loads the llh for the season.

                data = self.injectors[season["Name"]].create_dataset(0.0)
                llh = self.llhs[season["Name"]]

                # Loops over each source in catalogue

                for source in np.sorted(self.sources, order="Distance (Mpc)"):

                    # Identify spatially- and temporally-coincident data

                    mask = llh.select_spatially_coincident_data(data, [source])
                    spatial_coincident_data = data[mask]



                    t_mask = np.logical_and(
                        np.greater(spatial_coincident_data["timeMJD"],
                                   llh.time_pdf.sig_t0(source)),
                        np.less(spatial_coincident_data["timeMJD"],
                                llh.time_pdf.sig_t1(source))
                    )

                    coincident_data = spatial_coincident_data[t_mask]
                    total[j].append(len(coincident_data))
                    # If there are events in the window...

                    if len(coincident_data) > 0:

                        # Identify significant events (S/B > 1)

                        significant = llh.find_significant_events(
                            coincident_data, source)

                        results[j].append(len(significant))
                    else:
                        results[j].append(0)

        for j, season in enumerate(sorted(self.seasons)):
            res = results[j]
            tot = total[j]

            print season["Name"],"Significant events", np.mean(res), \
                np.median(res), np.std(res)
            print season["Name"], "All events", np.mean(tot), np.median(tot), \
                np.std(tot)

            llh = self.llhs[season["Name"]]

            for source in self.sources:

                print "Livetime", llh.time_pdf.effective_injection_time(source)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Path for analysis pkl_file")
    cfg = parser.parse_args()

    with open(cfg.file) as f:
        mh_dict = Pickle.load(f)

    mh = MinimisationHandler(mh_dict)
    mh.iterate_run(mh_dict["scale"], n_steps=mh_dict["n_steps"],
                   n_trials=mh_dict["n_trials"])
