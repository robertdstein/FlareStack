"""Script to compare the sensitivity for each TDE catalogues as a function of
injected spectral index. Rather than the traditional flux at 1 GeV,
Sensitivities are given as the total integrated fluence across all sources,
and as the corresponding standard-candle-luminosity.
"""
from __future__ import print_function
from __future__ import division
import numpy as np
from flarestack.data.icecube import diffuse_8_year, diffuse_10_year
from flarestack.utils.catalogue_loader import load_catalogue
from flarestack.analyses.agn_cores.shared_agncores import \
    agn_subset_catalogue, complete_cats_north, complete_cats_north, agn_catalogue_name, agn_subset_catalogue_north
from flarestack.core.minimisation import MinimisationHandler

from flarestack.cluster import analyse, wait_for_cluster

import logging

import os
import psutil, resource  #to get memory usage info

analyses = dict()

# Initialise Injectors/LLHs

llh_time = {
    "time_pdf_name": "Steady"
}

llh_energy = {
    "energy_pdf_name": "PowerLaw"
}

llh_dict = {
    "llh_name": "standard_matrix",
    "llh_sig_time_pdf": llh_time,
    "llh_energy_pdf": llh_energy
}


def base_name(cat_key, gamma):
    return "analyses/agn_cores/test_completeness_diffuse_8years_15steps_3scale_cluster_scaleby3_newcat/{0}/" \
           "{1}/".format(cat_key, gamma)


def generate_name(cat_key, n_sources, gamma):
    return base_name(cat_key, gamma) + "NrSrcs={0}/".format(n_sources)


gammas = [2.0, 2.5]

nr_brightest_sources = [1, 3, 10, 31, 100, 316, 1000, 3162, 10000, 32249]

all_res = dict()

running_time = []
for (cat_type, method) in complete_cats_north[1:2]:

    unique_key = cat_type + "_" + method

    print(unique_key)

    gamma_dict = dict()

    for gamma_index in gammas:
        res = dict()
        for j, nr_srcs in enumerate(nr_brightest_sources):

            cat_path = agn_subset_catalogue(cat_type, method, nr_srcs)
            catalogue = load_catalogue(cat_path)
            cat = np.load(cat_path)
            full_name = generate_name(unique_key, nr_srcs, gamma_index)

            injection_time = llh_time
            injection_energy = dict(llh_energy)
            injection_energy["gamma"] = gamma_index

            inj_kwargs = {
                "injection_energy_pdf": injection_energy,
                "injection_sig_time_pdf": injection_time,
            }

            mh_dict = {
                "name": full_name,
                "mh_name": "large_catalogue",
                "dataset": diffuse_8_year.get_seasons(), #subselection_fraction=1),
                "catalogue": cat_path,
                "llh_dict": llh_dict,
                "inj_dict": inj_kwargs,
                "n_trials": 1, #10,
            }

            mh = MinimisationHandler.create(mh_dict)

            # Set the scale for the injection
            scale_factor = 3 * mh.guess_scale()/3

            '''
            UNCOMMENT THIS IF:
            1. It is the first time you are running this code
            2. You want to run locally
            3. If you are running on the cluster with < 1000 sources
            '''
            mh_dict["n_steps"] = 15
            mh_dict["scale"] = scale_factor
            analyse(mh_dict, cluster=True, n_cpu=8, n_jobs=100,
                    len_catalogue=nr_srcs, nr_seasons=10)

            '''
            UNCOMMENT THIS IF:
            1. If you are running on the cluster with > 1000 sources
            '''
            # _n_jobs = 100
            # scale_loop = np.linspace(0, scale_factor, 15)
            # print(scale_loop)
            # for scale in scale_loop[:]:
            #     print('Running ' + str(mh_dict["n_trials"]) + ' trials with scale ' + str(scale))
            #     mh_dict["fixed_scale"] = scale
            #     if scale == 0.:
            #         n_jobs = _n_jobs*10
            #     else:
            #         n_jobs = _n_jobs
            #     print("Submitting " + str(n_jobs) + " jobs")
            #     analyse(mh_dict, cluster=True, n_cpu=1, n_jobs=n_jobs,
            #             len_catalogue=nr_srcs, nr_seasons=10)

            res[nr_srcs] = mh_dict

        gamma_dict[gamma_index] = res
    all_res[unique_key] = gamma_dict


wait_for_cluster()

logging.getLogger().setLevel("INFO")


for (cat_key, gamma_dict) in all_res.items():
    # print (cat_key, gamma_dict)
    agn_type = cat_key.split("_")[0]
    print(agn_type)
    xray_cat = cat_key.split(str(agn_type)+'_')[-1]
    print(xray_cat)

    full_cat = load_catalogue(agn_catalogue_name(agn_type, xray_cat))

    full_flux = np.sum(full_cat["base_weight"])

    saturate_ratio = 0.26

    for (gamma_index, gamma_res) in (iter(gamma_dict.items())):

        print("gamma: ", gamma_index)
        sens = []
        sens_err_low = []
        sens_err_upp = []
        disc_pot = []
        disc_ts_threshold = []
        sens_livetime = []
        n_src = []
        fracs = []
        disc_pots_livetime = []
        ratio_sens = []
        ratio_disc = []
        int_xray_flux_erg = []
        int_xray_flux = []
        guess = []
        sens_n = []
        disc_pot_n = []

        base_dir = base_name(cat_key, gamma_index)

        for (nr_srcs, rh_dict) in sorted(gamma_res.items()):

            cat = load_catalogue(rh_dict["catalogue"])

            print("nr_srcs in loop: ", nr_srcs)
            print("   ")
            print("   ")
            int_xray = np.sum(cat["base_weight"] / 1e13*624.151)
            int_xray_flux.append(int_xray) # GeV cm-2 s-1
            int_xray_flux_erg.append(np.sum(cat["base_weight"]) / 1e13) # erg
            # cm-2 s-1
            fracs.append(np.sum(cat["base_weight"])/full_flux)

            try:
                rh = ResultsHandler(rh_dict)
                print("Sens", rh.sensitivity)
                print("Sens_err", rh.sensitivity_err, rh.sensitivity_err[0], rh.sensitivity_err[1])
                print("Disc", rh.disc_potential)
                print("Disc_TS_threshold", rh.disc_ts_threshold)
                # print("Guess", rh_dict["scale"])
                print("Sens (n)", rh.sensitivity * rh.flux_to_ns)
                print("DP (n)", rh.disc_potential * rh.flux_to_ns)
                # # guess.append(k_to_flux(rh_dict["scale"])* 2./3.)
                # guess.append(k_to_flux(rh_dict["scale"])/3.)

                # sensitivity/dp normalized per flux normalization GeV-1 cm-2 s-1
                sens.append(rh.sensitivity)
                sens_err_low.append(rh.sensitivity_err[0])
                sens_err_upp.append(rh.sensitivity_err[1])
                disc_pot.append(rh.disc_potential)
                disc_ts_threshold.append(rh.disc_ts_threshold)
                sens_n.append(rh.sensitivity * rh.flux_to_ns)
                disc_pot_n.append(rh.disc_potential * rh.flux_to_ns)

                astro_sens, astro_disc = rh.astro_values(
                    rh_dict["inj_dict"]["injection_energy_pdf"])

                key = "Energy Flux (GeV cm^{-2} s^{-1})" # old version: "Total Fluence (GeV cm^{-2} s^{-1})"

                sens_livetime.append(astro_sens[key])   # fluence=integrated over energy
                disc_pots_livetime.append(astro_disc[key])

                ratio_sens.append(astro_sens[key] / int_xray) # fluence
                # normalized over tot xray flux
                ratio_disc.append(astro_disc[key] / int_xray)

                n_src.append(nr_srcs)

            except OSError:
                pass

        # n_src_modified = [1,3,10,31,100,1000,13927]
        # Save arrays to file
        np.savetxt(plot_output_dir(base_dir) + "data.out",
                   (np.array(n_src), np.array(int_xray_flux_erg), np.array(int_xray_flux),
                    np.array(sens), np.array(sens_err_low), np.array(sens_err_upp),
                    np.array(disc_pot), np.array(disc_ts_threshold),
                    np.array(sens_livetime), np.array(disc_pots_livetime),
                    np.array(ratio_sens), np.array(ratio_disc),
                    np.array(ratio_sens)/saturate_ratio, np.array(ratio_disc)/saturate_ratio,
                    np.array(sens_n), np.array(disc_pot_n)),
                    header="n_src, int_xray_flux_erg, int_xray_flux_gev," 
                           "sensitivity, sensitivity_err_lower, sensitivity_err_upper,"
                           "dp, disc_ts_threshold,"
                           "int_sensitivity, int_dp," 
                           "ratio_sens, ratio_dp,"
                           "ratio_sens_saturate, ratio_dp_saturate,"
                           "sensitivity_nr_neutrinos, dp_nr_neutrinos")

        labels = ['Sensitivity', 'Discovery Potential', 'sens', 'dp']
        for i, sens_dp in enumerate([sens_livetime, disc_pots_livetime]):
            print(i, sens_dp)
            # Plot 1: sensitivity/dp fluence vs number of sources (cuts)
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(n_src, sens_dp, ls='-', marker='o', label = labels[i])
            ax1.grid(True, which='both')
            ax1.set_ylabel(r"Total $nu$ flux [GeV cm$^{-2}$ s$^{-1}$]", fontsize=12)
            ax1.set_xlabel(r"Number of sources", fontsize=12)
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax1.set_xscale("log")
            ax1.set_yscale("log")
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] + "_vs_Nsrcs.pdf")
            plt.close()

        labels = ['Sensitivity', 'Discovery Potential', 'sens', 'dp']
        for i, sens_dp in enumerate([sens_livetime, disc_pots_livetime]):
            print(i, sens_dp)
            # Plot 2: sensitivity/dp fluence vs fraction of flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(fracs, sens_dp, ls='-', marker='o', label = labels[i])
            ax1.grid(True, which='both')
            ax1.set_ylabel(r"Total $nu$ flux [GeV cm$^{-2}$ s$^{-1}$]", fontsize=12)
            ax1.set_xlabel(r"Fraction of total X-Ray flux", fontsize=12)
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax1.set_xscale("log")
            ax1.set_yscale("log")
            ax1.set_xlim(1e-2, 1)
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] + "_vs_frac.pdf")
            plt.close()

        for i, sens_dp in enumerate([sens_livetime, disc_pots_livetime]):
            # Plot 3: sensitivity fluence vs integrated x-ray flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(int_xray_flux_erg, sens_dp, marker = 'o', ls='-',  label=labels[i])
            ax1.grid(True, which='both')
            # ax1.legend(loc='upper left', framealpha=1)
            ax1.set_ylabel(r"Total Fluence [GeV cm$^{-2}$ s$^{-1}$]", fontsize=12)
            ax1.set_xlabel(r"Integrated X-Ray Flux [erg cm$^{-2}$ s$^{-1}$]",  fontsize=12)
            ax1.set_xscale("log")
            ax1.set_yscale("log")
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] + "_vs_Xray_cut.pdf")
            plt.close()

        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):
            print(i, sens_dp)
            # Plot 1: sensitivity/dp fluence vs number of sources (cuts)
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(n_src, sens_dp, marker ='o', ls='', label ='Ratio ' + labels[i])
            ax1.grid(True, which='both')
            # ax1.legend(loc='upper left', framealpha=1)
            ax1.set_ylabel(r"Integrated $\nu$ / X-Ray flux ratio", fontsize=12)
            ax1.set_xlabel(r"Number of brightest sources",  fontsize=12)
            ax1.set_xscale("log")
            # ax1.set_yscale("log")
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] +
                        "_vs_Nsrcs_ratio.pdf")
            plt.close()

        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):
            print(i, sens_dp)
            # Plot 2: sensitivity/dp fluence vs fraction of total flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(fracs, sens_dp, marker ='o', ls='',
                     label='Ratio ' + labels[i])
            ax1.grid(True, which='both')
            # ax1.legend(loc='upper left', framealpha=1)
            ax1.set_ylabel(r"Integrated $\nu$ / X-Ray flux ratio", fontsize=12)
            ax1.set_xlabel(r"Fraction of total X-Ray flux", fontsize=12)
            ax1.set_xscale("log")
            # ax1.set_yscale("log")
            ax1.set_xlim(1e-2, 1)
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] +
                        "_vs_frac_ratio.pdf")
            plt.close()

        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):
            # Plot 3: sensitivity fluence vs integrated x-ray flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(int_xray_flux_erg, sens_dp, ls='', marker='o', label='Ratio ' + labels[i])
            ax1.grid(True, which='both')
            # ax1.legend(loc='upper left', framealpha=1)
            ax1.set_ylabel(r"Integrated $\nu$ / X-Ray flux ratio", fontsize=20)
            # ax1.set_ylabel(r"Total Fluence [GeV cm$^{-2}$ s$^{-1}$]", fontsize=12)
            ax1.set_xlabel(r"Integrated X-Ray Flux [erg cm$^{-2}$ s$^{-1}$]",  fontsize=12)
            ax1.set_xscale("log")
            # ax1.set_yscale("log")
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] +
                        "_vs_Xray_cut_ratio.pdf")
            plt.close()

        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):
            print(i, sens_dp)
            # Plot 2: sensitivity/dp fluence vs fraction of total flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(fracs, np.array(sens_dp) / saturate_ratio, marker ='o', ls='',
                     label='Ratio ' + labels[i])
            ax1.grid(True, which='both')
            # ax1.legend(loc='upper left', framealpha=1)
            ax1.set_ylabel(labels[i] + r" (cumulative population $\nu$ flux/diffuse $\nu$ flux)", fontsize=10)
            ax1.set_xlabel(r"Fraction of total X-Ray flux", fontsize=10)
            ax1.set_xscale("log")
            # ax1.set_yscale("log")
            ax1.set_ylim(0,3)
            ax1.set_xlim(1e-2, 1)
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            import matplotlib.ticker as ticker
            ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))

            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] +
                        "_vs_frac_diffuse.pdf")
            plt.close()

        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):
            # Plot 3: sensitivity fluence vs integrated x-ray flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)
            cols = ["#00A6EB", "#F79646", "g", "r"]
            ax1.plot(int_xray_flux_erg, np.array(sens_dp) / saturate_ratio, ls='', marker='o', label='Ratio ' + labels[i])
            ax1.grid(True, which='both')
            # ax1.legend(loc='upper left', framealpha=1)
            ax1.set_ylabel(labels[i] + r" (cumulative population $\nu$ flux/diffuse $\nu$ flux)", fontsize=12)
            ax1.set_xlabel(r"Integrated X-Ray Flux [erg cm$^{-2}$ s$^{-1}$]",  fontsize=12)
            ax1.set_xscale("log")
            # ax1.set_yscale("log")
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.tight_layout()
            plt.savefig(plot_output_dir(base_dir) + labels[i+2] +
                        "_vs_Xray_cut_diffuse.pdf")
            plt.close()


        # labels = ['Number of Sources', 'Integrated X-Ray Flux [erg cm$^{-2}$ s$^{-1}$]', 'Nsrcs', 'Xray_cut_ratio']
        # for i, xaxis in enumerate([n_src, int_xray_flux_erg]):
        #     print(i, sens_dp)
        #     # Plot 1: sensitivity/dp fluence vs number of sources (cuts)
        #     plt.rc('axes', axisbelow=True)
        #     plt.figure()
        #     ax1 = plt.subplot(111)
        #     ax1.plot(xaxis, guess, ls='-', marker='o', label = labels[i])
        #     ax1.grid(True, which='both')
        #     ax1.set_ylabel(r"Estimated Discovery Potential [GeV sr$^{-1}$ cm$^{-2}$ s$^{-1}$]", fontsize=12)
        #     ax1.set_xlabel(labels[i], fontsize=12)
        #     ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        #     ax1.set_xscale("log")
        #     ax1.set_yscale("log")
        #     plt.tight_layout()
        #     plt.savefig(plot_output_dir(base_dir) + "guess_vs_" +labels[i+2] +".pdf")
        #     plt.close()
        #
        # labels = ['Number of Sources', 'Integrated X-Ray Flux [erg cm$^{-2}$ s$^{-1}$]', 'Nsrcs', 'Xray_cut_ratio']
        # for i, xaxis in enumerate([n_src, int_xray_flux_erg]):
        #     print(i, sens_dp)
        #     # Plot 1: sensitivity/dp fluence vs number of sources (cuts)
        #     plt.rc('axes', axisbelow=True)
        #     plt.figure()
        #     ax1 = plt.subplot(111)
        #     ax1.plot(xaxis, np.array(guess)/np.array(disc_pots_livetime), ls='-', marker='o', label = labels[i])
        #     ax1.grid(True, which='both')
        #     ax1.set_ylabel(r"Guess / Discovery Potential ", fontsize=12)
        #     ax1.set_xlabel(labels[i], fontsize=12)
        #     ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        #     ax1.set_xscale("log")
        #     ax1.set_yscale("log")
        #     plt.tight_layout()
        #     plt.savefig(plot_output_dir(base_dir) + "ratio_guess_dp_vs_" +labels[i+2] +".pdf")
        #     plt.close()
        #
        # labels = ['Number of Sources', 'Integrated X-Ray Flux [erg cm$^{-2}$ s$^{-1}$]', 'Nsrcs', 'Xray_cut_ratio']
        # for i, xaxis in enumerate([n_src, int_xray_flux_erg]):
        #     print(i, sens_dp)
        #     # Plot 1: sensitivity/dp fluence vs number of sources (cuts)
        #     plt.rc('axes', axisbelow=True)
        #     plt.figure()
        #     ax1 = plt.subplot(111)
        #     ax1.plot(xaxis, np.array(guess)/np.array(sens_livetime), ls='-', marker='o', label = labels[i])
        #     ax1.grid(True, which='both')
        #     ax1.set_ylabel(r"Guess / Sensitivity", fontsize=12)
        #     ax1.set_xlabel(labels[i], fontsize=12)
        #     ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        #     ax1.set_xscale("log")
        #     ax1.set_yscale("log")
        #     plt.tight_layout()
        #     plt.savefig(plot_output_dir(base_dir) + "ratio_guess_sen_vs_" +labels[i+2] +".pdf")
        #     plt.close()


        labels = ['Sensitivity', 'Discovery Potential', 'sens', 'dp']
        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):

            sens_dp = np.array(sens_dp)
            int_xray_flux_erg = np.array(int_xray_flux_erg)
            n_src = np.array(n_src)

            mask_nan = sens_dp >= 0    # mask nan values
            int_xray_flux_erg_mask = int_xray_flux_erg[mask_nan]
            sens_dp = sens_dp[mask_nan]
            n_src_mask = n_src[mask_nan]

            # Plot 3: sensitivity fluence vs integrated x-ray flux
            plt.rc('axes', axisbelow=True)
            plt.figure()
            ax1 = plt.subplot(111)

            print("Making plot of DP/sensitivity")
            print (n_src, int_xray_flux_erg, np.array(sens_dp) / saturate_ratio)

            fig, ax1 = plt.subplots()
            ax1.plot(int_xray_flux_erg_mask, np.array(sens_dp) / saturate_ratio, ls='-', color='C0', lw=2,
                     # marker = 'o',
                     label= labels[i])
            ax1.set_xlabel(r"Integrated X-Ray flux of all stacked sources [erg cm$^{-2}$ s$^{-1}$]")
            ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax1.grid(True, which='both')
            # Make the y-axis label, ticks and tick labels match the line color.
            ax1.set_ylabel(r'Fraction of $\nu$ diffuse flux', color='C0')
            # ax1.tick_params('y', colors='C0')
            ax1.set_xscale("log")
            ax1.axhline(y=1, ls='--', lw=2, color='black', alpha=0.6)

            ax2 = ax1.twinx()
            ax2.plot(int_xray_flux_erg_mask, n_src_mask, ls='-', color='C1', lw=2,
                     # marker = 'o',
                     label='Number of stacked sources')
            ax2.set_ylabel('Number of stacked sources', color='C1')
            ax2.set_yscale("log")
            ax2.grid(False)
            # ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            # ax2.tick_params('y', colors='C1')

            plt.annotate("IceCube Preliminary ", (0.1, 0.9), alpha=0.7, fontsize=20, xycoords="axes fraction",
                         multialignment="center")
            fig.tight_layout()
            fig.savefig(plot_output_dir(base_dir) + labels[i + 2] +  "_vs_Xray_cut_diffuse_double_yaxis.pdf")
            plt.close('all')

        labels = ['Sensitivity', 'Discovery Potential', 'sens', 'dp']
        for i, sens_dp in enumerate([ratio_sens, ratio_disc]):
            sens_dp = np.array(sens_dp)
            int_xray_flux_erg = np.array(int_xray_flux_erg)
            n_src = np.array(n_src)
            mask_nan = sens_dp>=0    # mask nan values
            int_xray_flux_erg_mask = int_xray_flux_erg[mask_nan]
            sens_dp = sens_dp[mask_nan]
            n_src_mask = n_src[mask_nan]

            #PLOT
            # ax1 = plt.subplot(111)
            # fig1 = plt.figure(1, figsize=(10, 6.1803398875))
            fig1 = plt.figure(1, figsize=(10*1.4, 10*1.4/1.61))
            frame1 = fig1.add_axes((.1, .3, .8, .6))
            plt.plot(int_xray_flux_erg_mask, np.array(sens_dp) / saturate_ratio, ls='-', lw=2,)
                     # marker = 'o',)
            #          label= labels[i])
            frame1.set_xscale(u"log")
            frame1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            frame1.grid(True, which='both')
            plt.ylabel(r'Fraction of $\nu$ diffuse flux')
            plt.axhline(y=1, ls='--', lw=2, color='black', alpha=0.5)

            plt.annotate("IceCube Preliminary ", (0.1, 0.9), alpha=0.7, fontsize=22, xycoords="axes fraction",
                         multialignment="center")

            frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
            # frame1.tick_params(axis='y',which="major", direction='in', length=10)
            frame1.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off

            frame1.tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                left=False,  # ticks along the bottom edge are off
                right=False,)  # ticks along the top edge are off
                # labelbottom=False)  # labels along the bottom edge are off

            #Residual plot
            frame2=fig1.add_axes((.1,.1,.8,.19))
            plt.plot(int_xray_flux_erg_mask, n_src_mask, ls='-',  lw=2, color='C1',
                     # marker='o',
                     label='Number of stacked sources')
            frame2.set_xlabel(r"Integrated X-Ray flux of all stacked sources [erg cm$^{-2}$ s$^{-1}$]")
            # frame2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            frame2.set_xscale(u"log")
            frame2.set_yscale(u"log")
            frame2.set_ylabel('Stacked\nsources')
            frame2.set_yticks([10, 1000, 10000])

            from matplotlib.ticker import AutoMinorLocator
            from matplotlib.ticker import StrMethodFormatter
            from matplotlib import ticker
            minorLocator = AutoMinorLocator()
            # frame2.xaxis.set_minor_locator(minorLocator)
            # frame2.tick_params(axis='y', which='minor',direction = 'in',length = 5)
            # frame2.tick_params(which='major',direction = 'in',length = 10)
            frame2.grid(True, which='both')
            frame2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
            frame2.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,)  # ticks along the top edge are off
                # labelbottom=False)  # labels along the bottom edge are off

            frame2.tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                left=False,  # ticks along the bottom edge are off
                right=False,)  # ticks along the top edge are off
                # labelbottom=False)  # labels along the bottom edge are off

            fig1.tight_layout()
            fig1.savefig(plot_output_dir(base_dir) + labels[i + 2] + "_vs_Xray_cut_diffuse_double_plot.pdf")
            plt.close()



        # Sensitivity + DP
        sens_plot = np.array(ratio_sens)
        dp_plot = np.array(ratio_disc)

        mask_dp = dp_plot>=0    # mask nan values
        mask_sen = sens_plot >= 0

        int_xray_flux_erg_mask_dp = int_xray_flux_erg[mask_dp]
        int_xray_flux_erg_mask_sen = int_xray_flux_erg[mask_sen]

        sens_plot = sens_plot[mask_sen]
        dp_plot = dp_plot[mask_dp]
        n_src_mask_sens = n_src[mask_sen]

        #PLOT
        # fig1 = plt.figure(1, figsize=(10, 6.1803398875))
        fig1 = plt.figure(1, figsize=(10*1.4, 10*1.4/1.61))
        frame1 = fig1.add_axes((.1, .3, .8, .6))
        plt.plot(int_xray_flux_erg_mask_dp, dp_plot / saturate_ratio, ls='-', lw=2,
                 # marker = 'o',
                 label = r'Discovery Potential ($5\sigma$)')
        plt.plot(int_xray_flux_erg_mask_sen, sens_plot / saturate_ratio, color='C0',
                 ls='-.', lw=2,
                 # marker = 'o',
                 label= 'Sensitivity (90%)')

        frame1.set_xscale(u"log")
        # frame1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        frame1.grid(True, which='both')
        frame1.legend(loc='upper right', prop={'size':20}, framealpha=1)
        plt.ylabel(r'Fraction of $\nu$ diffuse flux')
        plt.axhline(y=1, ls='--', lw=2, color='black', alpha=0.5)

        plt.annotate("IceCube Preliminary ", (0.1, 0.9), alpha=0.7, fontsize=22, xycoords="axes fraction",
                     multialignment="center")
        # frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
        frame1.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

        frame1.tick_params(
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            left=False,  # ticks along the bottom edge are off
            right=False, )  # ticks along the top edge are off
        # labelbottom=False)  # labels along the bottom edge are off

        frame1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        ########################
        # Nr sources plot

        frame2=fig1.add_axes((.1,.1,.8,.19))
        plt.plot(int_xray_flux_erg_mask_sen, n_src_mask_sens, ls='-',  lw=2, color='C1',
                 # marker='o',
                 label='Number of stacked sources')
        frame2.set_xlabel(r"Integrated X-Ray flux of all stacked sources [erg cm$^{-2}$ s$^{-1}$]")
        # frame2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        frame2.set_xscale(u"log")
        frame2.set_yscale(u"log")
        frame2.set_yticks([10, 1000, 10000])
        frame2.set_ylabel('Stacked\nsources')

        from matplotlib.ticker import AutoMinorLocator
        from matplotlib.ticker import StrMethodFormatter
        from matplotlib import ticker
        minorLocator = AutoMinorLocator()
        # frame2.xaxis.set_minor_locator(minorLocator)
        # frame2.tick_params(axis='y', which='minor',direction = 'in',length = 5)
        # frame2.tick_params(which='major',direction = 'in',length = 10)
        frame2.grid(True, which='both')
        frame2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
        frame2.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,)   # ticks along the top edge are off
            # labelbottom=False)  # labels along the bottom edge are off

        frame2.tick_params(
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            left=False,  # ticks along the bottom edge are off
            right=False, )  # ticks along the top edge are off
        # labelbottom=False)  # labels along the bottom edge are off

        fig1.tight_layout()
        fig1.savefig(plot_output_dir(base_dir) + "sens_dp_vs_Xray_cut_diffuse_double_plot.pdf")
        plt.close()