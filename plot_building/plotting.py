import argparse
import ROOT
import numpy as np
import sys

# Import our own plotting functions
import plotting_functions as pf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", type=str,
                        help="Which rat production version would you like to run over? [Defaults to 5.0.1]",
                        default="5.0.1")
    parser.add_argument("-p", "--phase", type=str,
                        help="Which phase data should be considered? [Defaults to Solar]",
                        default="Solar")
    parser.add_argument("-d", "--data", type=str,
                        help="Path to parent data repository. [Defaults to /lustre/scratch/epp/neutrino/snoplus/production/hdf5_spectra_v0.3]",
                        default="/lustre/scratch/epp/neutrino/snoplus/production/hdf5_spectra_v0.3")
    parser.add_argument("-y", "--years", type=float,
                        help="Number of live years [Defaults to 0.5]",
                        default=0.5)
    args = parser.parse_args()

    # Define data path using input args
    data_path = "%s/%s/" % (args.data, args.version)

    # Get scalings
    internal, external, solar  =  pf.scalings_dict(liveYears=args.years)
    #total_scalings = dict( internal.items() + external.items() + solar.items())
    total_scalings = dict( internal.items() + solar.items())

    # Find spectra associated with each background
    #in_mc_paths, in_reco_paths, in_missing = pf.find_background_paths(internal, data_path, phase=args.phase)
    #ex_mc_paths, ex_reco_paths, ex_missing  = pf.find_background_paths(external, data_path, phase=args.phase)
    #solar_mc_paths, solar_reco_paths, solar_missing  = pf.find_background_paths(solar, data_path, phase=args.phase)

    # Concatenate interesting parameters 
    #mc_paths = dict( in_mc_paths.items() + ex_mc_paths.items() + solar_mc_paths.items())
    #reco_paths = dict( in_reco_paths.items() + ex_reco_paths.items() + solar_reco_paths.items())
    #scalings = dict( internal.items() + external.items() + solar.items())
    
    mc_paths, reco_paths, missing = pf.find_background_paths(total_scalings, data_path, phase=args.phase)

    # Create canvas
    can = ROOT.TCanvas('can', 'can')
    can.SetLogx(False)

    # Get stacks
    mc_sig_stack, mc_leg = pf.plot_signals(mc_paths, total_scalings, dimension='energy_mc')
    reco_sig_stack, reco_leg = pf.plot_signals(reco_paths, total_scalings, dimension='energy_reco')

    # Get sig and bkg sums
    mc_sig_sum = mc_sig_stack.GetStack().Last()
    reco_sig_sum = reco_sig_stack.GetStack().Last()

    # Draw MC plot & add summed lines
    mc_sig_sum.SetLineColor(1)
    mc_sig_sum.SetLineStyle(2)
    mc_sig_stack.Draw('nostack')
    mc_sig_sum.Draw('same')
    mc_sig_stack.SetTitle('Truth background spectrum')
    mc_sig_stack.GetXaxis().SetTitle('Energy [MeV]')
    mc_sig_stack.GetXaxis().SetLimits(0.2, 7)
    mc_sig_stack.GetYaxis().SetRangeUser(1., 1.e8)
    mc_leg.AddEntry(mc_sig_sum, 'Sum sig')
    mc_leg.Draw()
    can.SetLogx()
    can.SetLogy()
    can.Update()
    can.SaveAs('Backgrounds_mc.pdf')

    # Draw MC plot & add summed lines
    reco_sig_sum.SetLineColor(1)
    reco_sig_sum.SetLineStyle(2)
    reco_sig_stack.Draw('nostack')
    reco_sig_sum.Draw('same')
    reco_sig_stack.SetTitle('Reconstructed background spectrum')
    reco_sig_stack.GetXaxis().SetTitle('Energy [MeV]')
    reco_sig_stack.GetXaxis().SetLimits(0.2, 7)
    reco_sig_stack.GetYaxis().SetLimits(1., 1.e8)
    reco_leg.AddEntry(mc_sig_sum, 'Sum sig')
    reco_leg.Draw()
    can.SetLogx()
    can.SetLogy()
    can.Update()
    can.SaveAs('Backgrounds_reco.pdf')
