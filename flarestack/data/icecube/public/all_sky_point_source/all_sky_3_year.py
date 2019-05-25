import os
import numpy as np
import csv
import pickle
from flarestack.shared import public_dataset_dir, \
    energy_proxy_path, med_ang_res_path, effective_area_plot_path,\
    ang_res_plot_path
from flarestack.utils.make_SoB_splines import make_individual_spline_set
from flarestack.shared import SoB_spline_path, dataset_plot_dir
from flarestack.data import Dataset
from flarestack.data.icecube.public import PublicICSeason
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import zipfile

src_dir = os.path.dirname(os.path.realpath(__file__)) + "/"

zip_file = src_dir + "raw_data/3year-data-release.zip"

output_base_dir = public_dataset_dir + "all_sky_3_year/"
extract_dir = output_base_dir + "extracted_data"
data_dir = extract_dir + "/3year-data-release/"
output_data_dir = output_base_dir + "events/"
pseudo_mc_dir = output_data_dir + "pseudo_mc/"

# If data has not been extracted, then extract from zip file

if not os.path.isdir(data_dir):

    with zipfile.ZipFile(zip_file, "r") as zip_ref:
        zip_ref.extractall(extract_dir)


for path in [output_data_dir, pseudo_mc_dir]:

    try:
        os.makedirs(path)
    except OSError:
        pass


def data_path(season):
    return output_data_dir + season + ".npy"


def pseudo_mc_path(season):
    return pseudo_mc_dir + season + ".npy"

def pseudo_mc_binning(season):
    return pseudo_mc_dir + season + "_binning.npy"


data_dtype = np.dtype([
    ('ra', np.float),
    ('dec', np.float),
    ('logE', np.float),
    ('sigma', np.float),
    ('time', np.float),
    ('sinDec', np.float)
])


datasets = ["IC79-2010", "IC86-2011", "IC86-2012"]


def parse_numpy_dataset():
    """Function to parse the .txt file  of events into a numpy format
    readable by flarestack, which is the saved in the products/ subdirectory.
    """

    for dataset in datasets:

        data = []

        path = data_dir + dataset + "-events.txt"

        with open(path, "r") as f:

            csv_reader = csv.reader(f, delimiter=" ")

            for i, row in enumerate(csv_reader):
                if i > 0:
                    row = [float(x) for x in row if x != ""]

                    entry = (np.deg2rad(row[3]), np.deg2rad(row[4]),
                             row[1], np.deg2rad(row[2]),
                             row[0], np.sin(np.deg2rad(row[4]))
                             )

                    data.append(entry)

        data = np.array(data, dtype=data_dtype)

        exp_path = data_path(dataset)

        with open(exp_path, "wb") as f:
            print("Saving converted numpy array to", exp_path)
            pickle.dump(data, f)


sample_name = "all_sky_3_year"

ps_3_year = Dataset()


def make_season(season_name):
    season = PublicICSeason(
        season_name=season_name,
        sample_name=sample_name,
        exp_path=data_path(season_name),
        pseudo_mc_path=pseudo_mc_path(season_name),
        sin_dec_bins=np.linspace(-1., 1., 50),
        log_e_bins=np.arange(2., 9. + 0.01, 0.25)
    )
    ps_3_year.add_season(season)


for season_name in datasets:
    make_season(season_name)


# if __name__=="__main__":
#     parse_numpy_dataset()

def parse_angular_resolution():
    """Function to parse angular resolution."""
    for dataset in datasets:

        path = data_dir + dataset + "-AngRes.txt"

        x = []
        y = []

        with open(path, "r") as f:

            csv_reader = csv.reader(f, delimiter=" ")

            for i, row in enumerate(csv_reader):

                if i > 0:
                    row = [float(x) for x in row if x != ""]

                    true_e = 0.5*(row[0] + row[1])
                    log_e = np.log10(true_e)
                    med_ang_err = np.deg2rad(row[2])
                    x.append(log_e)
                    y.append(med_ang_err)

        x = np.array(x)
        y = np.array(y)

        # Kinematic angle accounts for ~ 1 degree / sqrt(E/TeV)
        # This is the angle between the neutrino and the reconstructed muon
        # This serves as a lower limit on neutrino angular resolution

        def kinematic_angle(log_e):
            e = 10**log_e
            return 1.0 * np.sqrt(10**3/e)

        z = np.linspace(1, 6, 100)

        plt.figure()
        plt.scatter(x, np.degrees(y), label="Published Median (Energy Proxy)",
                    color="orange")

        plt.plot(z, kinematic_angle(z),
                 label=r"$\nu-\mu$ Kinematic Angle (True Energy)")

        new_x = np.linspace(1, 3, 7)
        plt.scatter(new_x, kinematic_angle(new_x), color="purple",
                    label="Flarestack K.A. anchor values")

        # Remove Kinematic Angle region and recalculate angular resolution
        # using finer step size

        mask = x > 3.

        x = x[mask]
        y = y[mask]

        full_x = list(new_x) + list(x)
        full_y = list(np.deg2rad(kinematic_angle(new_x))) + list(y)

        plt.plot(full_x, np.degrees(full_y), color="red",
                 linestyle=":", label="Flarestack interpolation")

        plt.xlabel(r"$log_{10}(Energy)$")
        plt.ylabel("Median Angular Resolution (degrees)")
        plt.ylim(ymin=0)
        plt.legend()

        save_path = ang_res_plot_path(ps_3_year.seasons[dataset])

        try:
            os.makedirs(os.path.dirname(save_path))
        except OSError:
            pass

        plt.savefig(save_path)
        plt.close()

        ar_path = med_ang_res_path(ps_3_year.seasons[dataset])

        try:
            os.makedirs(os.path.dirname(ar_path))
        except OSError:
            pass

        with open(ar_path, "wb") as f:
            print("Saving converted numpy array to", ar_path)
            pickle.dump([full_x, full_y], f)


def parse_effective_areas():
    """Function to parse effective areas .txt into a format that flarestack
    can use to build Signal/Background splines.
    """
    data_dtype = np.dtype([
        ('logE', np.float),
        ('trueE', np.float),
        ('sinDec', np.float),
        ('trueDec', np.float),
        ('ow', np.float),
        ('a_eff', np.float),
        ("sigma", np.float)
    ])

    for dataset in datasets:

        pseudo_mc = []

        path = data_dir + dataset + "-TabulatedAeff.txt"

        exp = data_loader(data_path(dataset))
        lower_e = min(exp["logE"])
        upper_e = max(exp["logE"])

        # Select only upgoing muons. For these events, the dominant
        # background is atmospheric neutrinos with a known spectrum of E^-3.7.
        # Downgoing events, on the other hand, are contaminated by sneaking
        # muon bundles which are harder to model.

        exp = exp[exp["sinDec"] > 0.]

        with open(path, "r") as f:

            csv_reader = csv.reader(f, delimiter=" ")

            for i, row in enumerate(csv_reader):

                if i > 0:
                    row = [float(x) for x in row if x != ""]

                    true_e = 0.5*(row[0] + row[1])
                    log_e = np.log10(true_e)
                    sin_dec = -0.5*(row[2] + row[3])
                    true_dec = np.arcsin(sin_dec)
                    a_eff = row[4]

                    randoms = [log_e]

                    for log_e in randoms:

                        # if log_e < 3.:
                        #     factor = 1e-4
                        # else:
                        #     factor = 1.
                        factor = 1.

                        entry = tuple([
                            log_e, true_e, sin_dec, true_dec,
                            a_eff, a_eff, np.nan
                        ])

                        pseudo_mc.append(entry)

        pseudo_mc = np.array(pseudo_mc, dtype=data_dtype)

        plt.figure()
        ax1 = plt.subplot(311)
        res = ax1.hist(exp["logE"], density=True)

        exp_vals = res[0]
        exp_bins = res[1]
        ax1.set_yscale("log")
        ax2 = plt.subplot(312, sharex=ax1)
        res = ax2.hist(
            pseudo_mc["logE"],
            weights=pseudo_mc["ow"] * pseudo_mc["trueE"]**-3.7,
            density=True, bins=exp_bins)
        mc_vals = res[0]

        ax2.set_yscale("log")

        # Maps ratio of expected neutrino energies to energy proxy values
        # This can tell us about how true energy maps to energy proxy

        centers = 0.5 * (exp_bins[:-1] + exp_bins[1:])

        # Fill in empty bins

        mc_vals = np.array(mc_vals)

        # print(mc_vals)

        mc_vals += min(pseudo_mc["ow"][pseudo_mc["ow"] > 0.]) * centers ** -3.7
        # print(mc_vals)

        x = [-5.0] + list(centers) + [15.0]
        y = exp_vals / mc_vals
        y = [y[0]] + list(y) + [y[-1]]

        log_e_weighting = interp1d(x, np.log(y))

        ax3 = plt.subplot(313)
        plt.plot(centers, exp_vals/mc_vals)
        plt.plot(centers, np.exp(log_e_weighting(centers)),
                 linestyle=":")
        ax3.set_yscale("log")

        save_path = effective_area_plot_path(ps_3_year.seasons[dataset])

        try:
            os.makedirs(os.path.dirname(save_path))
        except OSError:
            pass

        plt.savefig(save_path)
        plt.close()

        pseudo_mc["ow"] *= np.exp(log_e_weighting(pseudo_mc["logE"]))

        mc_path = pseudo_mc_path(dataset)

        with open(mc_path, "wb") as f:
            print("Saving converted numpy array to", mc_path)
            pickle.dump(pseudo_mc, f)

        ep_path = energy_proxy_path(ps_3_year.seasons[dataset])

        try:
            os.makedirs(os.path.dirname(ep_path))
        except OSError:
            pass

        with open(ep_path, "wb") as f:
            print("Saving converted numpy array to", ep_path)
            pickle.dump([x, np.log(y)], f)


if __name__ == "__main__":
    from flarestack.icecube_utils.dataset_loader import data_loader

    # mc = data_loader(ps_7year[0]["mc_path"])
    # print(mc.dtype.names)
    # for x in mc:
    #     true_e = x["trueE"]
    #     print(true_e, np.log10(true_e), x["logE"])
    #     input("prompt")
    parse_angular_resolution()
    parse_effective_areas()

    for season in ps_3_year.get_seasons().values():
        # exp = data_loader(season["exp_path"])
        # mc = data_loader(season["pseudo_mc"])
        # make_individual_spline_set(season, SoB_spline_path(season))
        season.plot_effective_area()
