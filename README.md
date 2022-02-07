# get_sea_spray_mode.m
 A MATLAB function to retrieve the sea spray aerosol mode from measured submicron size distributions and supermicron scattering.
 
This code applies a Mie inversion to retrieve lognormal fitting parameters of the sea spray mode using submicron size distributions and supermicron scattering. The code requires the sea spray Mie scattering look-up table (**sea_spray_mie_table.mat**) to be in the same directory as the retrieval code (or accessible from a working directory).
 
 
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
