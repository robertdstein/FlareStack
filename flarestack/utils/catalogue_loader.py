import numpy as np
from numpy.lib.recfunctions import append_fields, rename_fields


def calculate_source_weight(sources) -> float:
    """
    Depending on the dimension of the input, calculate:
    - the sum of the weights for a vector of sources
    - the absolute weight of an individual source.

    Dividing the absolute weight of a single source by the sum of the weights of all sources in the catalogue gives the relative weight of the source, i.e. the fraction of the fitted flux that is produced by it. The fraction of the injected flux may be different, since the "injection weight modifier" is not accounted for.
    """
    return np.sum(sources["base_weight"] * sources["distance_mpc"] ** -2)


def get_relative_source_weight(catalogue, source):
    """
    Get the relative weight of a source.
    """
    return calculate_source_weight(source) / calculate_source_weight(catalogue)


def load_catalogue(path):

    sources = np.load(path)

    # Maintain backwards-compatibility

    maps = [
        ("ra", "ra_rad"),
        ("dec", "dec_rad"),
        ("Relative Injection Weight", "injection_weight_modifier"),
        ("Distance (Mpc)", "distance_mpc"),
        ("Name", "source_name"),
        ("Ref Time (MJD)", "ref_time_mjd"),
        ("Start Time (MJD)", "start_time_mjd"),
        ("End Time (MJD)", "end_time_mjd"),
    ]

    for (old_key, new_key) in maps:

        if old_key in sources.dtype.names:

            sources = rename_fields(sources, {old_key: new_key})

    if "base_weight" not in sources.dtype.names:

        base_weight = np.ones(len(sources))

        sources = append_fields(
            sources, "base_weight", base_weight, usemask=False, dtypes=[float]
        )

    # Check that ra and dec are really in radians!

    if max(sources["ra_rad"]) > 2.0 * np.pi:
        raise Exception(
            "Sources have Right Ascension values greater than 2 "
            "pi. Are you sure you're not using degrees rather "
            "than radians?"
        )

    if max(abs(sources["dec_rad"])) > np.pi / 2.0:
        raise Exception(
            "Sources have Declination values exceeding "
            "+/- pi/2. Are you sure you're not using degrees "
            "rather than radians?"
        )

    # Check that all sources have a unique name

    if len(set(sources["source_name"])) < len(sources["source_name"]):

        raise Exception(
            "Some sources in catalogue do not have unique "
            "names. Please assign unique names to each source."
        )

    # Rescale 'base_weight'
    # sources["base_weight"] /= np.mean(sources["base_weight"])

    # Order sources
    sources = np.sort(sources, order="distance_mpc")

    return sources


# def convert_catalogue(path):
#     print "Converting", path
#     cat = load_catalogue(path)
#     print cat
#     # np.save(path, cat)

# if __name__ == "__main__":
#     import os
#     from flarestack.analyses.agn_cores.shared_agncores import agncores_cat_dir
#
#     # for path in os.listdir(agncores_cat_dir):
#     for path in ["radioloud_radioselected_100brightest_srcs.npy"]:
#         filename = agncores_cat_dir + path
#         cat = load_catalogue(filename)
#         cat["base_weight"] = cat['injection_weight_modifier']
#         cat['injection_weight_modifier'] = np.ones_like(cat["base_weight"])
#         np.save(filename, cat)
