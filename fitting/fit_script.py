'''An example script for fitting systematics for a single 
background signal. All other Backgrounds have been held
as fixed

Author: Ed leming <e.leming@sussex.ac.uk>
Date  : 03/03/2016
'''
import argparse
import json 
import time
import os
import sys
import ROOT

# Local imports
import plot_building.plotting_functions as pf

# Echidna imports
import echidna.core.spectra as spectra
import echidna.core.config as config
import echidna.core.smear as smear
import echidna.core.scale as scale
import echidna.core.shift as shift
import echidna.output.store as store
import echidna.output.plot_root as plot_root
import echidna.fit.test_statistic as test_statistic
import echidna.fit.fit as fit


def generate_spec_fit_config(prior, spectra_name, sigma=0.05, low=0, high=0, bins=11):
    '''Generate a SpectrumFitConfig object for signal
    '''
    sigma = sigma*prior
    if low == 0:
        low = prior - (10*sigma)
    if high == 0:
        high = prior + (10*sigma)
    rate_dict = {}
    rate_dict["prior"] = prior
    rate_dict["sigma"] = sigma
    rate_dict["low"] = low
    rate_dict["high"] = high
    rate_dict["bins"] = bins
    spec_dict = {}
    spec_dict["rate"] = rate_dict
    config_dict = {}
    config_dict["spectral_fit_parameters"] = spec_dict
    '''
    config_dict = ({"spectral_fit_parameters": { 
                       "rate": {
                          "prior": prior,
                          "sigma": sigma, 
                          "low"  : low,
                          "high" : high,
                          "bins" : bins}}})
    '''
    config_object = config.SpectraFitConfig.load(config_dict, spectra_name, name="%s_SpectFitConfig" % spectra_name)
    config_object.dump_to_file("tmp_")
    config_object = config.SpectraFitConfig.load_from_file("tmp_spectral_fit_config.yml", spectra_name)
    return config_object


def smear_mc(spectrum):
    ''' Smear the mc spectrum to SNO+'s expected energy resolution
    '''
    smearer = smear.EnergySmearLY()
    smearer.set_num_sigma(3)
    #smearer.set_resolution(0.05)
    return smearer.weighted_smear(spectrum)
    

def shrink(spectrum, roi):
    ''' Shrink the spectrum to a specified roi
    '''
    # Do it this way so works for both mc and reco
    pars =  spectrum.get_config().get_pars()
    shrink_dict = {"%s_low" % pars[0]: roi["energy"][0],
                    "%s_high" % pars[0]: roi["energy"][1],
                    "%s_low" % pars[1]: roi["radial3"][0],
                    "%s_high" % pars[1]: roi["radial3"][1]}
    spectrum.shrink(**shrink_dict)
    return spectrum

def make_plot(data_spec, fixed_spec, floated_spec, name, initial=True):
    ''' Make a root plot to show what's going on
    '''
    can = ROOT.TCanvas("c1","c1")
    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    sig_stack = ROOT.THStack()

    data = plot_root.plot_projection(data_spec, "energy_reco", graphical=False)
    fixed_backgrounds = plot_root.plot_projection(fixed_spec, "energy_reco", graphical=False)
    floated_background = plot_root.plot_projection(floated_spec, "energy_mc", graphical=False)

    # Characteristic stuff
    data.SetLineColor(1)
    fixed_backgrounds.SetLineColor(12)
    floated_background.SetLineColor(2)
    # Make backgrounds dotted
    fixed_backgrounds.SetLineStyle(2)

    # Make legend 
    leg.AddEntry(data, "data")
    leg.AddEntry(fixed_backgrounds, "fixed_backgrounds")
    leg.AddEntry(floated_background, name)

    # Add to the stack
    sig_stack.Add(data)
    sig_stack.Add(fixed_backgrounds)
    sig_stack.Add(floated_background)

    sig_stack.Draw('nostack')
    if initial:
        sig_stack.SetTitle('Initial conditions for fitting %s' % name)
    else:
        sig_stack.SetTitle('Best fit result for %s' % name)
    sig_stack.GetXaxis().SetTitle('Energy [MeV]')
    sig_stack.GetXaxis().SetLimits(0, 5)
    #sig_stack.GetYaxis().SetRangeUser(1., 1.e8)
    leg.Draw()
    can.SetLogy()
    can.Update()
    raw_input("Hit enter to return")
    if initial:
        can.SaveAs('results/initial_conditions/%s.pdf' % name)
    else:
        can.SaveAs('results/best_fit/%s.pdf' % name)
    # Tidy up a little
    del can
    del sig_stack

def apply_best_fit(spectrum, fit_summary, name):
    ''' Apply the systematics associated with the best fit
    '''
    # Shift
    shifter = shift.Shift()
    shifter.set_shift(fit_summary["energy_mc_shift"]["best_fit"])
    spectrum = shifter.shift(spectrum, "energy_mc")

    # Scale
    scaler = scale.Scale()
    scaler.set_scale_factor(fit_summary["energy_mc_scale"]["best_fit"])
    spectrum = scaler.scale(spectrum, "energy_mc")
    
    # Apply rate - confusingly called scale
    spectrum.scale(fit_summary["%s_rate" % name]["best_fit"])
    return spectrum

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", type=str,
                        help="Which rat production version would you like to run over? [Defaults to 5.0.1]",
                        default="5.0.1")
    parser.add_argument("-p", "--phase", type=str,
                        help="Which phase data should be considered? [Defaults to Solar]",
                        default="Solar")
    parser.add_argument("-b", "--bkgnd", type=str,
                        help="Background (or `module') spectrum to be fitted for")
    parser.add_argument("-y", "--years", type=float,
                        help="Number of live years [Defaults to 0.5]",
                        default=0.5)
    args = parser.parse_args()
    
    # Check user defined background. If included, remove the phase
    if args.bkgnd[:len(args.phase)] == args.phase:
        signal = args.bkgnd[len(args.phase):]
    else:
        signal = args.bkgnd

    # load in background scalings 
    internal, external, solar = pf.scalings_dict()
    #scale_dict = dict( internal.items() + solar.items())
    scale_dict = dict( internal.items() + solar.items() + external.items())
    # Check passed background exists in dict    
    
    # Find paths to spectra files
    data_path = "/lustre/scratch/epp/neutrino/snoplus/production/hdf5_spectra_v0.3/%s/" % args.version
    mc_paths, reco_paths, missing = pf.find_background_paths(scale_dict, data_path, phase=args.phase)

    # Check if signal path exists before continuing
    try:
        mc_paths[signal]
    except KeyError:
        print "Spectrum file does not exist for: %s!" % signal
        sys.exit()

    # Load spec config from file
    this_dir = os.getcwd()
    spec_config = config.SpectraConfig.load_from_file("%s/data_config.yml" % this_dir)

    #################################################
    # Create data spectra and fixed backgrounds dict
    #################################################
    # Define roi for fit.
    roi = {"energy": (0.0, 6.0), "radial3": (0., 0.2)}
    # Create empty spectra and create binning vars
    n_energy_bins, n_rad_bins = int(250*(roi["energy"][1]/10.)), int(200*roi["radial3"][1])
    data = spectra.Spectra("data", 0, spec_config)
    data = shrink(data, roi)
    data.rebin((n_energy_bins, n_rad_bins))
    fixed_backgrounds = spectra.Spectra("fixed_background", 0, spec_config)
    fixed_backgrounds = shrink(fixed_backgrounds, roi)
    fixed_backgrounds.rebin((n_energy_bins, n_rad_bins))
    fixed_backgrounds_dict = {}
    for idx, key in enumerate(scale_dict.keys()):
        # Check key exists
        try:
            reco_paths[key]
        except KeyError:
            print 'Warning: Spectrum file does not exist for: %s ' % key
            continue

        # Load spectrum
        spec = store.load(reco_paths[key])
        # Shrink spectra to roi
        spec = shrink(spec, roi)
        spec.rebin((n_energy_bins, n_rad_bins))
        # Scale
        spec.scale(scale_dict[key])
        # Add reco spectrum to the 'data' spectra
        data.add(spec)

        # Add spectra to data and append fixed bkgnds
        if key == signal:
            # Loadd mc spec for fitting
            signal_spec = store.load(mc_paths[key])
            signal_spec = shrink(signal_spec, roi)
            signal_spec.rebin((n_energy_bins, n_rad_bins))
            signal_spec.scale(scale_dict[key])
            signal_spec = smear_mc(signal_spec)
            spec_fit_config = generate_spec_fit_config(scale_dict[key], signal, bins=11)
            signal_spec.set_fit_config(spec_fit_config)
            print "\nLoaded %s as signal\n" % key
        else:
            # Add to fixed backgrounds
            fixed_backgrounds.add(spec)
            fixed_backgrounds_dict[spec] = spec.get_num_decays()
            print "Loaded %s as fixed background" % key

    print "Spectra sucessfully loaded"
    make_plot(data, fixed_backgrounds, signal_spec, signal, initial=True)

    ###########################################
    # Create fitter
    ###########################################
    fit_config = config.GlobalFitConfig.load_from_file("%s/fit_config.yml" % this_dir)
    chi_squared = test_statistic.BakerCousinsChi(per_bin=True)

    # Create fitter
    fitter = fit.Fit(roi, chi_squared, fit_config, data=data,
                     fixed_backgrounds=fixed_backgrounds_dict,
                     floating_backgrounds=[signal_spec], per_bin=True, 
                     use_pre_made=False)
    
    # Run fitter
    print "Varying systematics.... This will take a few minutes"
    stat_zero = fitter.fit()
    fit_results = fitter.get_fit_results()
    ###########################################
    # Print and save results
    ###########################################
    print "RESULTS: "
    print "\n%s\n" % json.dumps(fit_results.get_summary())

    # Apply best fit
    final_spec = apply_best_fit(signal_spec, fit_results.get_summary(), signal)
    make_plot(data, fixed_backgrounds, final_spec, signal, initial=False)

    # Save
    timestr = time.strftime("%d_%m_%Y-%H.%M")
    resultsFile = "%s/results/%s_%s.dat" % (this_dir, signal, timestr)
    print "Saving data to: %s" % resultsFile
    store.dump_fit_results(resultsFile, fit_results)
