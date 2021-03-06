# Undergrad project for Echidna

A brief libaray of functions and scripts to run a fit in echidna. Paths to both mc and reco echidna spectrum objects
(saved as .hdf5) are hardcoded assuming the user is running on the Sussex 'feynman' cluster. 

## plot_building

plotting.py: A script to plot a background spectrum. It takes arguments: -v, -p, -d, and -y to define the production version
	     [5.0.1], phase [Solar], data path [feynman] and number of years of running [0.5]. This script loads spectra 
	     objects associated with the passed production version and phase and builds a background spectrum from the 'reco'
	     spectra. 

plotting_functions.py: Some functions called by plotting.py. Undergrads we asked to edit the functions to perform various tasks.

plot_fit_results.py: An example script to plot the results generated by an echidna fitting script. The undergrads would
		     use this as an example to generate plots relavent to their analysis.

## fitting

fit_script.py: A script to perform an example echina fit. It loads reco spectra from file and sums them into a 'signal' spectrum.
	       A single mc spectra is loaded for the isotope of interest (defined with the -b flag) and smeared to have an 
	       energy resolution, as given by poission statistics at 200 Ph/MeV. The energy shift, energy scale and rate 
	       systematics are floated and chi2 values taken at each step. The results of each chi2 test are stored in the
	       ./results directory, which will be autogenerated. Plots of the inital and best fit parameters for the smeared
	       'background' spectrum is also saved in the ./results direc.