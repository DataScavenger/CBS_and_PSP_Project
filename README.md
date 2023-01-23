This repository contains the scripts and functions used for the paper "Slowing of Frontal Beta Oscillations in Atypical Parkinsonism".


--------
Software
--------

Software used for the data analysis:

- Matlab R2018a
- Fieldtrip (ver. 14.12.2020)
- Python 3.9.1
- FOOOF - fitting oscillations & one over f (development ver. 1.0.1)
- colorbrewer2 (ver. 18.09.2015)
- ROInets (part of OSL: https://ohba-analysis.github.io/osl-docs/)
- RainCloudPlots (ver. 12.08.2018)
- drawbrace (https://de.mathworks.com/matlabcentral/fileexchange/38716-curly-brace-annotation)

-------
Scripts
-------

Run the following analysis scripts in succession. Selfcoded functions used for the script are indented.

__I Preprocessing
- cbs_sensor2areas_info
- cbs_prep_databrowse
		- cbs_cutstring
		- cbs_waterfall
- cbs_prep_ica_calculation
- cbs_prep_ica_visualize
- cbs_prep_cleaning
- cbs_prep_headmodel

__II Parcel-Analysis					
- cbs_source_parcellation_aal_cloud
	- cbs_sumareas
- cbs_parcel_visualization
- cbs_source_reconstruction
- cbs_parcel_time2power
- cbs_parcel_fooof_ana_single_flexible.py
- cbs_parcel_fooof_paper_example.py
- cbs_parcel_fooof_power
	- cbs_clean_subjects
	- cbs_mirror_activity_parcel
- cbs_parcel_fooof_power_statistics
	- cbs_prepare_neighbours
- cbs_parcel_fooof_CoM
- cbs_parcel_fooof_CoM_statistics
	- cbs_prepare_neighbours
- cbs_parcel_fooof_CoM_Crr
- cbs_parcel_fooof_spectral_features_group_I
- cbs_parcel_fooof_spectral_features_statistics
- cbs_source_power_individual
