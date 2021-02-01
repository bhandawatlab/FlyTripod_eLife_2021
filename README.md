# FlyTripod_eLife_2021
Codes by Chanwoo Chun <cc2465@cornell.edu>

Description for each subfolder in MATLAB folder:

1. Acquisition:	Codes for automatically acquiring video data during experiments.
		DrosoLoco_Auto.m is the main function

2. Tracking: 	Codes for post-processing the acquired video data.
		Tracks CoM and leg tips (saves to shot.mat data files)
		mainTrack.m is the main function.

3. Models:	makeStepData.m segmentalize all shot.mat data into individual tripod steps.
		For each tripod step (stance), fit*.m files fits mechanical models to the step data.
		Outputs stepData.m (stepData*.m) file

4. Plotting:	Codes for plotting figures
		The codes use shot.mat and stepData files.
		
		Figure Number:	Code Name:
		Fig 2D		durations.m
		Fig 2-S2 	fluctuation_trend.m
		Fig 3A		compositeGaitmap_final.m
		Fig 3B,C,E,F 	cycle_analysis.m
		Fig 4		fluctuation_trend.m
		Fig 6A,B,C	ARSLIP_Plots.m
		Fig 6-S1 	fitted_params.m
		Fig 7C	 	SpringyTripod.m
		Fig 6D; 7D,E 	ARSLIP_Plot_speedControl.m
		Fig 8 		SpringyTripod.m

5. Functions:	Contains custom functions


Description for each subfolder in Data folder:

1. RawVideoExample: Contains a set of example raw data files saved by DrosoLoco_Auto.m

2. ProcessedData: Contains shot.m files and metadata from all experiments.

3. StepData: Contains example stepData*.m files
