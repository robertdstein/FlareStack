import numpy as np
import os
import cPickle as Pickle
from core.minimisation import MinimisationHandler
from core.results import ResultsHandler
from data.icecube_pointsource_7_year import ps_7year
from shared import plot_output_dir, flux_to_k, analysis_dir, catalogue_dir
from utils.skylab_reference import skylab_7year_sensitivity
from cluster import run_desy_cluster as rd
from utils.custom_seasons import custom_dataset
import matplotlib.pyplot as plt

name = "analyses/tde/compare_flare_to_box/"

analyses = dict()

cat_path = catalogue_dir + "TDEs/individual_TDEs/Swift J1644+57_catalogue.npy"
# cat_path = catalogue_dir + "TDEs/individual_TDEs/XMMSL1 J0740-85_catalogue.npy"
catalogue = np.load(cat_path)

t_start = catalogue["Start Time (MJD)"]
t_end = catalogue["End Time (MJD)"]

max_window = float(t_end - t_start)

# Initialise Injectors/LLHs

injection_energy = {
    "Name": "Power Law",
    "Gamma": 2.0,
}

llh_time = {
    "Name": "FixedEndBox",
}

llh_energy = injection_energy

no_flare = {
    "LLH Energy PDF": llh_energy,
    "LLH Time PDF": llh_time,
    "Fit Gamma?": True,
    "Flare Search?": False,
    "Fit Negative n_s?": False
}

no_flare_negative = {
    "LLH Energy PDF": llh_energy,
    "LLH Time PDF": llh_time,
    "Fit Gamma?": True,
    "Flare Search?": False,
    "Fit Negative n_s?": True
}

flare = {
    "LLH Energy PDF": llh_energy,
    "LLH Time PDF": llh_time,
    "Fit Gamma?": True,
    "Flare Search?": True,
    "Fit Negative n_s?": False
}

src_res = dict()

lengths = np.logspace(-2, 0, 9) * max_window
# lengths = np.logspace(-2, 0, 17) * max_window


for i, llh_kwargs in enumerate([
                                no_flare,
                                no_flare_negative,
                                flare
                                ]):

    label = ["Time-Integrated", "Time-Integrated (negative n_s)",
             "Cluster Search"][i]
    f_name = ["fixed_box", "fixed_box_negative", "flare"][i]

    flare_name = name + f_name + "/"

    res = dict()

    for flare_length in lengths:

        full_name = flare_name + str(flare_length) + "/"

        injection_time = {
            "Name": "FixedRefBox",
            "Fixed Ref Time (MJD)": t_start,
            "Pre-Window": 0,
            "Post-Window": flare_length,
            "Time Smear?": True,
            "Min Offset": 0.,
            "Max Offset": max_window - flare_length
        }

        inj_kwargs = {
            "Injection Energy PDF": injection_energy,
            "Injection Time PDF": injection_time,
            "Poisson Smear?": True,
        }

        scale = flux_to_k(skylab_7year_sensitivity(np.sin(catalogue["dec"]))
                          * (50 * max_window/ flare_length))

        mh_dict = {
            "name": full_name,
            "datasets": custom_dataset(ps_7year, catalogue,
                                       llh_kwargs["LLH Time PDF"]),
            "catalogue": cat_path,
            "inj kwargs": inj_kwargs,
            "llh kwargs": llh_kwargs,
            "scale": scale,
            "n_trials": 1,
            "n_steps": 15
        }

        analysis_path = analysis_dir + full_name

        try:
            os.makedirs(analysis_path)
        except OSError:
            pass

        pkl_file = analysis_path + "dict.pkl"

        with open(pkl_file, "wb") as f:
            Pickle.dump(mh_dict, f)

        # rd.submit_to_cluster(pkl_file, n_jobs=12000)

        # mh = MinimisationHandler(mh_dict)
        # mh.iterate_run(mh_dict["scale"], mh_dict["n_steps"], n_trials=3)
        # mh.clear()
        res[flare_length] = mh_dict

    src_res[label] = res

# rd.wait_for_cluster()

sens = [[] for _ in src_res]
fracs = [[] for _ in src_res]
disc_pots = [[] for _ in src_res]
sens_e = [[] for _ in src_res]
disc_e = [[] for _ in src_res]

labels = []

for i, (f_type, res) in enumerate(sorted(src_res.iteritems())):

    if f_type!="Time-Integrated (negative n_s)":
        for (length, rh_dict) in sorted(res.iteritems()):
            try:
                rh = ResultsHandler(rh_dict["name"], rh_dict["llh kwargs"],
                                    rh_dict["catalogue"], show_inj=True)

                inj_time = length * (60 * 60 * 24)

                astro_sens, astro_disc = rh.astro_values(
                    rh_dict["inj kwargs"]["Injection Energy PDF"])

                key = "Total Fluence (GeV cm^{-2} s^{-1})"

                e_key = "Mean Luminosity (erg/s)"

                sens[i].append(astro_sens[key] * inj_time)
                disc_pots[i].append(astro_disc[key] * inj_time)

                sens_e[i].append(astro_sens[e_key] * inj_time)
                disc_e[i].append(astro_disc[e_key] * inj_time)

                fracs[i].append(length)

            except OSError:
                pass

            except KeyError:
                pass

            except EOFError:
                pass

        labels.append(f_type)

for j, [fluence, energy] in enumerate([[sens, sens_e],
                                      [disc_pots, disc_e]]):

    plt.figure()
    ax1 = plt.subplot(111)

    ax2 = ax1.twinx()

    cols = ["#00A6EB", "#F79646", "g", "r"]
    linestyle = ["-", "-"][j]

    for i, f in enumerate(fracs):

        if len(f) > 0:

            ax1.plot(f, fluence[i], label=labels[i], linestyle=linestyle,
                     color=cols[i])
            ax2.plot(f, energy[i], linestyle=linestyle,
                     color=cols[i])

    ax2.grid(True, which='both')
    ax1.set_ylabel(r"Total Fluence [GeV cm$^{-2}$]", fontsize=12)
    ax2.set_ylabel(r"Mean Isotropic-Equivalent $E_{\nu}$ (erg)")
    ax1.set_xlabel(r"Flare Length (Days)")
    ax1.set_yscale("log")
    ax2.set_yscale("log")

    for k, ax in enumerate([ax1, ax2]):
        y = [fluence, energy][k]

        ax.set_ylim(0.95 * min([min(x) for x in y if len(x) > 0]),
                    1.1 * max([max(x) for x in y if len(x) > 0]))

    plt.title(["Sensitivity", "Discovery Potential"][j] + " for " + \
              catalogue["Name"][0])

    ax1.legend(loc='upper left', fancybox=True, framealpha=0.)
    plt.tight_layout()
    plt.savefig(plot_output_dir(name) + "/flare_vs_box_" +
                ["sens", "disc"][j] + ".pdf")
    plt.close()
