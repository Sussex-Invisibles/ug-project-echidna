import argparse
import numpy as np
import ROOT
import json

import echidna.output.store as store
from echidna.fit.fit_results import FitResults

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str,
                        help="Results .txt file for plotting")
    args = parser.parse_args()

    # Load fit results from file
    fit_results = store.load_fit_results(args.file)

    # Print summary
    print "%s" % json.dumps(fit_results.get_summary())
    
    fit_config = fit_results._fit_config
    spec_config = fit_results._spectra_config
    fit_pars = fit_results._fit_config.get_pars()
    spec_pars = fit_results._spectra_config.get_pars()

    # Rate par is named relatively to background.
    # Find it in the list of pars to make projection
    # easy.
    for par in fit_pars:
        if "rate" in par:
            rate_par = par

    # Get systematics values
    scale_values = fit_config.get_par("energy_mc_scale").get_values()
    shift_values = fit_config.get_par("energy_mc_shift").get_values()
    rate_values = fit_config.get_par(rate_par).get_values()
    energy_values = spec_config.get_par("energy_mc").get_bin_centres()
    radius_values = spec_config.get_par("radial3_mc").get_bin_centres()

    # Get chi2 projected onto a single axis
    scale = fit_results.nd_project_stats("energy_mc_scale")
    shift = fit_results.nd_project_stats("energy_mc_shift")
    rate = fit_results.nd_project_stats(rate_par)
    energy  = fit_results.nd_project_stats("energy_mc")
    radius = fit_results.nd_project_stats("radial3_mc")

    can = ROOT.TCanvas("c1", "c1")
    graph = ROOT.TGraph(len(scale), scale_values, scale)
    graph.SetTitle("Chi2 projection for scale values")
    graph.Draw("AP*")
    can.Update()
    raw_input("Enter to return")

    graph = ROOT.TGraph(len(shift), shift_values, shift)
    graph.SetTitle("Chi2 projection for shift values")
    graph.Draw("AP*")
    can.Update()
    raw_input("Enter to return")

    graph = ROOT.TGraph(len(rate), rate_values, rate)
    graph.SetTitle("Chi2 projection for rate values")
    graph.Draw("AP*")
    can.Update()
    raw_input("Enter to return")

    graph = ROOT.TGraph(len(energy), energy_values, energy)
    graph.SetTitle("Chi2 projection for energy values")
    graph.Draw("AP*")
    can.Update()
    raw_input("Enter to return")

    graph = ROOT.TGraph(len(radius), radius_values, radius)
    graph.SetTitle("Chi2 projection for radius values")
    graph.Draw("AP*")
    can.Update()
    raw_input("Enter to return")
