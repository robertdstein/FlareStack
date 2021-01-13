"""simple script to test the scale estimation implemented in the Submitter class"""

import logging
import unittest
from flarestack.shared import flux_to_k
from flarestack.data.public import icecube_ps_3_year
from flarestack.utils.prepare_catalogue import ps_catalogue_name
from flarestack.cluster.submitter import Submitter, DESYSubmitter

injection_energy = {
    "energy_pdf_name": "power_law",
    "gamma": 2.0,
}

injection_time = {
    "time_pdf_name": "steady",
}

llh_time = {
    "time_pdf_name": "steady",
}

inj_kwargs = {
    "injection_energy_pdf": injection_energy,
    "injection_sig_time_pdf": injection_time,
}

llh_energy = injection_energy

llh_kwargs = {
    "llh_name": "standard",
    "llh_energy_pdf": llh_energy,
    "llh_sig_time_pdf": llh_time,
    "llh_bkg_time_pdf": {"time_pdf_name": "steady"}
}

base_name = "test/test_submitter/"

sindec = 0.
cat_path = ps_catalogue_name(sindec)
scale = 0.39370132 * 5
mh_dict = {
    "name": base_name,
    "mh_name": "fixed_weights",
    "dataset": icecube_ps_3_year,
    "catalogue": cat_path,
    "inj_dict": inj_kwargs,
    "llh_dict": llh_kwargs,
    "scale": scale,
    "n_steps": 3
}

public_sens_3yr = 4.533328532314386e-10
upper = 7.66510624e-12
lower = 7.93338706e-12


class TestSubmitter(unittest.TestCase):

    def setUp(self):
        pass

    def test_submitter(self):
        logging.info('testing Submitter class')
        this_mh_dict = dict(mh_dict)
        this_mh_dict['name'] += 'test_submitter/'
        this_mh_dict['n_trials'] = 10
        sb = Submitter.get_submitter(this_mh_dict, use_cluster=False, n_cpu=5, remove_old_results=True)
        sb.analyse()

    def test_scale_estimation(self):
        this_mh_dict = dict(mh_dict)
        this_mh_dict['name'] += 'test_scale_estimation/'
        this_mh_dict['scale'] *= 5.1
        this_mh_dict['n_steps'] = 6
        sb = Submitter.get_submitter(this_mh_dict, use_cluster=False, n_cpu=5,
                                     do_sensitivity_scale_estimation='quick_injections')
        sb.run_quick_injections_to_estimate_sensitivity_scale()
        true_value = flux_to_k(public_sens_3yr)
        self.assertAlmostEqual(sb.sens_guess, true_value / 0.9, delta=true_value / 0.9 * 0.6)
        self.assertGreater(sb.sens_guess / 0.5, true_value)

    def test_desy_submitter(self):
        this_mh_dict = dict(mh_dict)
        this_mh_dict['name'] += 'test_desy_submitter/'
        desy_sb = DESYSubmitter(this_mh_dict, use_cluster=False, n_cpu=5)
        desy_sb.make_cluster_submission_script()


if __name__ == '__main__':
    unittest.main()