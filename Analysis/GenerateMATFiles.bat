:: Convert a list of DLC tracking data into MATLAB format. This must be run with the conda labeleditor environment active.

:: Set the conda environment
set CONDAPATH=C:\ProgramData\Anaconda2
set ENVNAME=labeleditor
if %ENVNAME%==base (set ENVPATH=%CONDAPATH%) else (set ENVPATH=%CONDAPATH%\envs\%ENVNAME%)

:: Activate the environment
call %CONDAPATH%\Scripts\activate.bat %ENVPATH%

@echo off
for %%x in (
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\165DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\185DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\203DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\208DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\215DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\216DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\230DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\241DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\245DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\248DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\250DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\259DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\260DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\263DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\266DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180321Sideview\274DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\64DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\65DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\66DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\67DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\72DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\75DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
		"D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\80DLC_resnet50_walkchamberJul7shuffle1_200000.h5"
       ) do (
         echo DLC file: %%x
		 python "D:\GitHub\dlc-label-editor\matconverter.py" %%x
         echo.
         echo =-=-=-=-=-=
         echo.
       )
	   
:: Deactivate the environment
:: call conda deactivate

pause
