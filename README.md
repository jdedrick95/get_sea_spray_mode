# get_sea_spray_mode.m
A MATLAB function to retrieve the sea spray aerosol mode from measured submicron size distributions and supermicron scattering.
 
This code applies a Mie inversion to retrieve lognormal fitting parameters of the sea spray mode using submicron size distributions and supermicron scattering. The code requires the sea spray Mie scattering look-up table (**sea_spray_mie_table.mat**) to be in the same directory as the retrieval code (or accessible from a working directory).

Methodology for this code is described in Dedrick et al. "Retrieval of the Sea Spray Aerosol Mode from Submicron Particle Size Distributions and Supermicron Scattering during LASIC" (submitted to *Atmospheric Measurement Techniques*)
 
**INPUTS:**

- **bsca_RGB**: supermicron scattering coefficients (Mm-1); row dimension: time; column dimension: [450nm, 550nm, 700nm]

- **bsca_std**: standard deviation of supermicron scattering coefficients (Mm-1) during temporal average; row dimension: time; column dimension: [450nm, 550nm, 700nm]

- **bsca_inst_std**: instrument scattering error/uncertainty (%), single value (e.g. 5% [Frie and Bahreini, 2021])

- **PNSD**: submicron particle size distribution (cm-3 µm-1); row dimension: concentration; column dimension: time

- **PNSD_D**: submicron particle diameters from size distribution (µm); row dimension: diameters

- **D_op**: overlap region for which to constrain Mie solutions (this variable is a user specified selection of diameters from PNSD_D) row dimension; diameters

- **PNSD_std**: standard deviation of particle size distribution during temporal average; resolved at each size bin (cm-3 µm-1); row dimension: concentration; column dimension: time 

- **PNSD_N_std**: Instrument concentration error/uncertainty (%), single value (e.g. 10% [Frie and Bahreini, 2021])

- **PNSD_D_std**: Instrument sizing uncertainty (%), single value (e.g. 2.5%, DMT UHSAS sizing uncertainty)

**OUTPUTS:**

- **sea_spray_mode**: sea spray mode fitting parameters; row dimension: time; column dimension: [number, mean diameter, geometric standard deviation] 

- **sea_spray_mode_95**: 95% confidence interval ranges of the sea spray mode fitting parameters. This variable is a cell matrix; row dimension: [number, mean diameter, geometric standard deviation]; columns dimension: [lower 95th, upper 95th]; cell dimension: time

- **error_thresh**: scattering error threshold (Mm-1); row dimension: time

- **low_error_idx**: indices of the look-up table that fall below the error threshold. This variable is a cell matrix; row dimension: Mie look-up table indices; cell dimension: time 

- **test_coeff**: probable Mie solutions that are tested against the measured size distribution. This variable is a cell matrix; row dimension: probable Mie solutions; cell dimension: time

- **RSS_fit**: residual sum of squares of unique sea spray mode to measured size distribution; row dimension: time

- **chi2_fit**: chi-square error of unique sea spray mode to measured size distribution; row dimension: time

- **D_mie**: size distribution diameters (µm)

- **dlogDp_mie**: log-base 10 difference of the diameters

- **fail_flag**: flag value identifying reason for retrieval failure (0 = retrieval successful, 1 = scattering not available at all 3 wavelengths, 2 = no Mie scattering solutions below the error threshold, 3 = no Mie solutions that are within the joint probability 95th percentile that can be tested against the size distribution); row dimension: time

- **retrieval_duration**: time to complete the retrieval (minutes); row dimension: time

# test.m

- test tex
