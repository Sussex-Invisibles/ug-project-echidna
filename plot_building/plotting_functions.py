""" Functions for plotting production data.
"""
import echidna.output.store as store
import echidna.output.plot_root as plot_root
import numpy as np
import ROOT
import argparse
import os
import sys

def scalings_dict(liveYears=0.5):
    '''Returns a dictionary containing background scalings associated with the number of live years
    passed. Background rates taken from Backgrounds_target_LAB.txt, associated with docDB 507v27

    Args:
        liveYears: the numer of liveYears the detector will have been running.
     
    Returns:
        scalingDict: a dictionary containing the standard scalings associated with each background isotope.
    '''
    ########################
    # Internal backgrounds
    ########################
    # Uranium chain
    internal = {}
    internal['U232'] = 4898
    internal['Th234'] = 4898
    internal['Pa234m'] = 4898
    internal['U234'] = 4898
    internal['Th230'] = 4898
    internal['Ra226'] = 4898
    internal['Rn222'] = 4898
    internal['Po218'] = 4898
    internal['Pb214'] = 4898
    internal['Bi212'] = 4898
    internal['Bi214'] = 4898
    internal['Po214'] = 4897
    internal['Tl210'] = 1
    internal['Pb210'] = 42671
    internal['Bi210'] = 42671
    internal['Po210'] = 1.7e7

    # Thorium chain
    internal['Th232'] = 682
    internal['Ra228'] = 682
    internal['Ac228'] = 682
    internal['Th228'] = 682
    internal['Ra224'] = 682
    internal['Rn220'] = 682
    internal['Po216'] = 682
    internal['Pb212'] = 682
    internal['Bi212'] = 682
    internal['Po212'] = 436
    internal['Tl208'] = 246

    #Others
    internal['Ar39'] = 85483
    internal['C14'] = 4.08e9
    internal['Kr85'] = 85483
    internal['K40'] = 8479
    # Scale for live years
    for key in internal:
        internal[key] = internal[key]*liveYears

    #########################################
    # Solar - Numbers from Richie Bonventre
    # bonventre@gmail.com
    #########################################
    solar = {}
    solar['B8_Solar'] = 1270
    solar['Be7_Solar'] = 139245
    solar['Cno_Solar'] = 31204
    solar['Pep_Solar'] = 7920
    for key in solar:
        solar[key] = solar[key]*liveYears

    ####################################
    # Externals.. A little more involved
    ####################################
    external = {}
    # Internal calib ropes
    external['Bi214_Intropes'] = 4966
    external['Tl210_Intropes'] = 1.04
    external['Tl208_Intropes'] = 418
    external['Bipo212-500_Intropes'] = 743
    external['K40_Intropes'] = 2.81e4

    # Hold down ropes
    external['Bi214_Hdropes'] = 4.06e6
    external['Tl210_Hdropes'] = 853
    external['Tl208_Hdropes'] = 2.32e6
    external['Bipo212-500_Hdropes'] = 4.13e6
    external['K40_Hdropes'] = 1.89e8

    # Hold up ropes
    external['Bi214_Hupropes'] = 8.35e5
    external['Tl210_Hupropes'] = 175
    external['Tl208_Hupropes'] = 4.78e5
    external['Bipo212-500_Hupropes'] = 8.5e5
    external['K40_Hupropes'] = 3.9e7

    # Water shielding -- !FOR RAT > 4.6!
    external['Bi214_Exwater'] = 1.32e8
    external['Tl210_Exwater'] = 2.77e4
    external['Tl208_Exwater'] = 3.92e6
    external['Bipo212-500_Exwater'] = 6.96e6

    # Acrylic vessel
    external['Bi214_Av'] = 1.28e7
    external['Tl210_Av'] = 2682
    external['Tl208_Av'] = 1.5e6
    external['Bipo212-500_Av'] = 2.67e6
    external['K40_Av'] = 7.32e7

    # AV dust, external
    external['Bi214_Av_Outerdust'] = 7.75e5
    external['Tl210_Av_Outerdust'] = 163
    external['Tl208_Av_Outerdust'] = 4.6e5
    external['Bipo212-500_Av_Outerdust'] = 8.2e5 
    external['K40_Av_Outerdust'] = 1.76e7

    # AV dust, internal
    external['Bi214_Av_Innerdust'] = 4.15e4
    external['Tl210_Av_Innderdust'] = 8.7
    external['Tl208_Av_Innerdust'] = 2.48e4
    external['Bipo212-500_Av_Innerdust'] = 4.41e4
    external['K40_Av_Innerdust'] = 9.4e5

    # PMTs
    external['Bi214_Pmt'] = 3.7e11
    external['Tl210_Pmt'] = 7.8e7
    external['Tl208_Pmt'] = 4.4e10
    external['Bipo212_Pmt'] = 7.8e10
    # Pmt_Betagammas???

    # AV radon
    external['Pb210_Av_Radon'] = 3.63e10
    external['Po210_Av_Radon'] = 3.63e10
    external['Bi210_Av_Radon'] = 3.63e10

    for key in external:
        external[key] = external[key]*liveYears

    ################################
    # Leaching stuff
    ################################
    Pb_av, Bi_av, Po_av, Pb_lab, Bi_lab, Po_lab = leaching_model(liveYears=liveYears)
    internal["Pb210"] = internal["Pb210"] + Pb_lab
    internal["Bi210"] = internal["Bi210"] + Bi_lab
    internal["Po210"] = internal["Po210"] + Po_lab
    print "%1.2e'\t%1.2e" % (Pb_av, Pb_lab)
    return internal, external, solar

def find_background_paths(scalings_dict, data_path, phase='TeLoaded'):
    '''Find hdf5 files associated with background keys defined in scalings_dict. Files are located 
    within the passed data_path
    '''
    data_dirs = os.listdir(data_path)
    mc_paths, reco_paths, missing = {}, {}, []
    for key in scalings_dict:
        for da in data_dirs:
            if da[len(phase):] ==  key:
                fList = os.listdir( '%s%s/Summed/' % (data_path, da) )
                for f in fList:
                    if 'reco_sum.hdf5' in f:
                        reco_paths[key] = '%s%s/Summed/%s' % (data_path, da, f)
                    elif 'mc_sum.hdf5' in f:
                        mc_paths[key] = '%s%s/Summed/%s' % (data_path, da, f)
        # If we couldn't find it, record the fact
        if not key in reco_paths:
            missing.append("%s%s" % (phase, key))
        elif not key in mc_paths:
            missing.append("%s%s" % (phase, key))
    return mc_paths, reco_paths, missing 

def leaching_model(Pb210_av=985, Bi210_av=985.62, Po210_av=1013.2, liveYears=0.5):
    '''Returns the number of counts at both the LAB and acrylic surface after some number of liveyears. The
    values for the initial activity of Pb210/Bi210/Po210 came from Valentina [Bq]. The generic leaching model was 
    written by Oleg. The basic formulars are:
    
    dPb_av = Pb210_av - Pb210_av_decayed - Pb210_leached
    dBi_av = Bi210_av - Bi210_av_decayed - Bi210_leached + Pb210_av_decayed
    dPo_av = Po210_av - Po210_av_decayed - Po210_leached + Bi210_av_decayed
    
    dPb_lab = Pb210_lab - Po210_lab_decayed + Pb210_leached
    dBi_lab = Bi210_lab - Bi210_lab_decayed + Bi210_leached + Pb210_lab_decayed
    dPo_lab = Po210_lab - Po210_lab_decayed + Po210_leached + Bi210_lab_decayed


    Returns:
        Pb210_lab [float]: The integrated number of counts leached into the lab in t=liveYears
        Bi210_lab [float]: The integrated number of counts leached into the lab in t=liveYears
        Po210_lab [float]: The integrated number of counts leached into the lab in t=liveYears
        Pb210_av  [float]: The integrated number of counts at the av surface in t=liveYears
        Bi210_av  [float]: The integrated number of counts at the av surface in t=liveYears
        Po210_av  [float]: The integrated number of counts at the av surface in t=liveYears
    '''
    # Convert liveYears into other units:
    days = int(np.ceil(liveYears*365.2))
    hours = int(np.ceil(liveYears*365.2*24))
    minutes = int(np.ceil(liveYears*365.2*24*60))

    one_day = 60*60*24 # [s]
    one_hour = 60*60   # [s]
    one_minute = 60    # [s]

    # Half lives
    half_life_Pb = 22.3*365.2*24*60*60 # [s] 
    half_life_Bi = 5.01*24*60*60       # [s] 
    half_life_Po = 138.367*24*60*60    # [s] 

    # Calculate decay rate (lambda) for different isotopes
    rate_Pb = np.log(2)/half_life_Pb   # [1/s]
    rate_Bi = np.log(2)/half_life_Bi   # [1/s]
    rate_Po = np.log(2)/half_life_Po   # [1/s]
    
    # Leaching rate: 1.4e-4 +/- 1.74e-5 day^-1 @ 12degC 
    # Leaching time: 7200 +/- 900 days 
    LR = 1.4e-4 / (60*60*24) #s^-1 

    # Convert activity to N_atoms
    atoms_Pb = Pb210_av / rate_Pb
    atoms_Bi = Bi210_av / rate_Bi
    atoms_Po = Po210_av / rate_Po

    ###############################################################
    # Calculate how many atoms of each isotope after some time, dt
    ###############################################################
    Pb_av, Bi_av, Po_av = np.zeros(minutes+1), np.zeros(minutes+1), np.zeros(minutes+1)
    Pb_lab, Bi_lab, Po_lab = np.zeros(minutes+1), np.zeros(minutes+1), np.zeros(minutes+1)
    # Set initial values at the av surface - The LAB activity is initially zero
    Pb_av[0], Bi_av[0], Po_av[0] = atoms_Pb, atoms_Bi, atoms_Po
    # Loop over the number of minutes we want to run for
    for i in range(minutes):
        # Calculate how many atoms decayed
        Pb210_av_decayed = Pb_av[i]*(1-np.exp(-rate_Pb*one_minute))
        Bi210_av_decayed = Bi_av[i]*(1-np.exp(-rate_Bi*one_minute))
        Po210_av_decayed = Po_av[i]*(1-np.exp(-rate_Po*one_minute))

        Pb210_lab_decayed = Pb_lab[i]*(1-np.exp(-rate_Pb*one_minute))
        Bi210_lab_decayed = Bi_lab[i]*(1-np.exp(-rate_Bi*one_minute))
        Po210_lab_decayed = Po_lab[i]*(1-np.exp(-rate_Po*one_minute))
        
        # Calculate how many atoms leached
        Pb210_leached = Pb_av[i]*(1-np.exp(-LR*one_minute))
        Bi210_leached = Bi_av[i]*(1-np.exp(-LR*one_minute))
        Po210_leached = Po_av[i]*(1-np.exp(-LR*one_minute))
        
        # Calculate how many atoms of each we'll have after some time, dt
        Pb_av[i+1] = Pb_av[i] - Pb210_av_decayed - Pb210_leached
        Bi_av[i+1] = Bi_av[i] - Bi210_av_decayed - Bi210_leached + Pb210_av_decayed
        Po_av[i+1] = Po_av[i] - Po210_av_decayed - Po210_leached + Bi210_av_decayed
        
        Pb_lab[i+1] = Pb_lab[i] - Pb210_lab_decayed + Pb210_leached
        Bi_lab[i+1] = Bi_lab[i] - Bi210_lab_decayed + Bi210_leached + Pb210_lab_decayed
        Po_lab[i+1] = Po_lab[i] - Po210_lab_decayed + Po210_leached + Bi210_lab_decayed
    
    # Convert N_atoms to Activity and integrate under
    # curve to give total counts.
    Pb_av_counts = np.trapz(Pb_av*rate_Pb, dx=one_minute)
    Bi_av_counts = np.trapz(Bi_av*rate_Bi, dx=one_minute)
    Po_av_counts = np.trapz(Po_av*rate_Po, dx=one_minute)

    Pb_lab_counts = np.trapz(Pb_lab*rate_Pb, dx=one_minute)
    Bi_lab_counts = np.trapz(Bi_lab*rate_Bi, dx=one_minute)
    Po_lab_counts = np.trapz(Po_lab*rate_Po, dx=one_minute)
    
    return Pb_av_counts, Bi_av_counts, Po_av_counts, Pb_lab_counts, Bi_lab_counts, Po_lab_counts


def plot_signals(spectra_dict, scale_dict,  dimension='energy_reco'):
    '''Pass a dictionary of paths to production hdf5 spectra and a associated dicitionary
    with scaling information. All signal spectra are then plotted individually, and a summed
    'signal' plotted as a thick black line.
    '''
    # Set-up stuff
    leg = ROOT.TLegend(0.5,0.6,0.9,0.9);
    leg.SetNColumns(2)
    sig_stack = ROOT.THStack()
    n_energy_bins, n_rad_bins = 250, 200
    # Loop over backgrounds and add to stack
    for idx, key in enumerate(scale_dict.keys()):
        try:
            spectra_dict[key]
        except KeyError:
            print 'Warning :Spectrum file does not exist for: %s ' % key
            continue
        # Load spectrum, scale and plot
        spec = store.load(spectra_dict[key])
        if key == 'B8_Solar':
            spec.shrink_to_roi(0, 10, dimension)
        # Rebin - need to make sure en and rad bins are correct way around.
        pars = spec.get_config().get_pars()
        if "energy" in pars[0]:
            spec.rebin([n_energy_bins, n_rad_bins])
        else:
            spec.rebin([n_rad_bins, n_energy_bins])
        spec.scale(scale_dict[key])
        plot = plot_root.plot_projection(spec, dimension, graphical=False)

        ######################################
        # ADD CUTS HERE
        cuts = run_cuts(plot)

        ######################################
        # DO FORMATTING HERE
        plot.SetLineColor(idx+2) # This can be replaced

        # plot is a TH1D object. So we can set line color etc.
        # 
        # Some thoughts....
        # 1) Are you sure you want to add this signal to the stack?
        #    We're interested in signals which we might be able to
        #    identify within the 'sum' response you see on the plot.
        #    Some signals will be so small we'll never be able to fit
        #    their signals with all those other backgrounds around.

        # 2) Does this signal belong to a group that we could perhaps 
        #    colour coordinate (Thorium chain, Uranium chain). It
        #    might make the final plot a little more read-able. 

        # 3) Are we missing any of the signals you picked out for
        #    you linearities? In an ideal world we would fit for those
        #    signals and try to reproduce the linearites you've 
        #    already produced

        # Update stack / legend
        if cuts:
            leg.AddEntry(plot, key)
            sig_stack.Add(plot)

    return sig_stack, leg

def run_cuts(plot):
    '''Run some basic cuts on the plot. Return True is this plot should be included.
    '''
    maxi = plot.GetMaximum()
    if maxi < 1:
        return False
    return True

def get_TColours(nColours):
    '''Defines new TColor objects  evenly distributed between red and blue with nColours steps
    Access them with ROOT.gROOT.GetColor(n) where n is between 1000 and 1000+nColours.
    '''
    cols = []
    step = 1. / nColours
    for i in range(nColours):
        r, g, b = rgb(0, 1, step*i)
        # Set to floats in range 0-1
        r = float(r / 255.)
        g = float(g / 255.)
        b = float(b / 255.)
        cols.append(ROOT.TColor(1000+i, r, g, b))
    return cols
    
def rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r
    return r, g, b    

def rescale_spectra_sum(spectrum, scaling):
    '''Rescale the *sum counts* of within a echida spectrum object. This is different
    to the scale() method of the spectra object which scales to the number of *decays*
    which would produce the spectrum.
    '''
    np.multiply(spectrum._data, scaling / spectrum.sum())
    return spectrum

def check_dir(dname): 
    '''Check if directory exists, if not then make it
    '''
    direc = os.path.dirname(dname)
    try:
        os.stat(direc)
    except:
        os.makedirs(direc)
        print "Made directory %s...." % dname

def main():
    '''Main is empty
    '''
    return 0

if __name__ == "__main__": 
    #Do nothing
    print main()
