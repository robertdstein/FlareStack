"""A standard time-integrated analysis is performed, using one year of
IceCube data (IC86_1).
"""
import logging
import unittest
import numpy as np
from flarestack.data.public import icecube_ps_3_year
from flarestack.utils.prepare_catalogue import ps_catalogue_name
from flarestack.core.unblinding import create_unblinder

# Initialise Injectors/LLHs

llh_dict = {
    "llh_name": "fixed_energy",
    "llh_sig_time_pdf": {
        "time_pdf_name": "steady"
    },
    "llh_energy_pdf": {
        "energy_pdf_name": "power_law",
        "gamma": 3.0
    }
}

# Loop over sin(dec) values

sindecs = np.linspace(0.5, -0.5, 3)


# These results arise from high-statistics sensitivity calculations,
# and can be considered the "true" answers. The results we obtain will be
# compared to these values.

true_parameters = [
    [0.0],
    [0.0],
    [2.0817298059169915]
]


class TestTimeIntegrated(unittest.TestCase):

    def setUp(self):
        pass

    def test_declination_sensitivity(self):
        logging.info("Testing 'fixed_energy' LLH class")

        # Test three declinations

        for j, sindec in enumerate(sindecs):

            unblind_dict = {
                "mh_name": "fixed_weights",
                "dataset": icecube_ps_3_year.get_seasons("IC86-2011"),
                "catalogue": ps_catalogue_name(sindec),
                "llh_dict": llh_dict,
            }

            ub = create_unblinder(unblind_dict)
            key = [x for x in ub.res_dict.keys() if x != "TS"][0]
            res = ub.res_dict[key]
            self.assertEqual(list(res["x"]), true_parameters[j])

            logging.info("Best fit values {0}".format(list(res)))
            logging.info("Reference best fit {0}".format(true_parameters[j]))


if __name__ == '__main__':
    unittest.main()
