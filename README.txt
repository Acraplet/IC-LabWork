This folder contains the code for the analysis of the data in the lab.

As of December 2024 this code analyses the optical power meter (OPM) data and Monitor PMT data.

OPM data (in the form of .csv file) should be stored in data/OPM/. The monitor PMT data (including the laser syncout data) should be stored in data/MonitorPMT.

## Analysis of OPM data:

Individual OPM files are analysed with the python/analyse_powermeter.py. Note that this code assumes the name of the file to extract a run number (I have been using data_powermeter_and_monitor_<laser_pulse_width>_<measurement_index>) which I have chosen as the laser pulse width for simplicity. This code should be called using the dedicated bash script: python/run_all_OPM_data_analysis.sh which calls the analyse_powermeter code for each of the files in a given folder.

This code then scrapes the .csv file (which have a standardised format) reading in the measurement date, the sampling rate, the wavelength measured by the OPM and the power read out at each sampling point.
It then plots a histogram of all of the entries as well as the optical power as a function of time. The Histogram is fitted with a gaussian, whose parameters are written out in the legend. These plots are saved in plots/OPM_per_run under (New_)OPM_meas_<laser_pulse_width>_<measurement_index>. It might be necessary to change this naming convention.

All the information about the run (including the fit results) are saved in a .text file in results/OMP_per_run and a global summary text file holding information about all of the processed OPM files is made under the name results/(new_)pulse_power_pd_summary.txt.

All of the informations stored in the summary .txt file is then read by the code python/plot_OPM_scans.py which makes the summary plots showing the power as a function of the pulse width and the PD current. Note that the conversion between the PD current and the pusle width is constant and not set by the user. The summary OPM plots are stored in /plots/summary.

## Analysis of Monitor PMT data:
The Monitor PMT data is analysed using the code/plot_volts_peaks executable that is compiled from the code/plot_volts_peaks.cpp code using the code/compile_analysis_code.sh script. Please remember to compile before running if you make any changes to the code. The current version of the code is reading in all of the files that are available by

The plot_volts_peaks executable reads in the root file which has two trees (in the case of the Monitor PMT only anlysis, there will be more when we look at the laser ball data), one called Trigger which corresponds to channel 2 on the oscilloscope and should always be connectected to the laser sync-out and one called PMT1Data holding the monitor PMT data (connected to Channel1 on the oscilloscope). The read-out voltage is stored for both trees in the branch called volts.

This code reads in each entry for the laser sync out (square wave) and uses that to identify rising edges to detetermine the start and end of each laser period. It then looks for the peak value in the monitor PMT within that period and stores it in the peaks histogram. It also sums the value of the data point before the peak, the peak value and the two points following the peak and stores that information in a separate histogram, this is the so called integrated peak. Finally it also makes a histogram with all of the recorded values (aka spectrum) of the recorded light into a third histogram. A gaussian is fitted to the distribution of integrated peaks and the fit mean and sigma, as well as the min and max values present in the waveform and the pulse width in a separate .txt file at results/Monitor_PMT_volts_summary.txt.

Please note, due to the high number of laser periods in the sample this process is pretty slow, we might consider speeding it up some other way for lower intensity laser ball runs if necessary.

This code outputs .pdf plots in the ../plots/summary/ folder with the spectrum, peak and integrated peaks distributions. It also saves a .root file under the name /results/Monitor_PMT_volts_overlay_spectrum.root with all of the histograms for each of the runs considered in case we need to perform additional anlaysis on it.


##Combination between OPM and Monitor PMT data:
After we have run the OPM and Monitor PMT data analysis separately we can combine them to validate the linearity of the monitor PMT. To do so we use the code python/plot_Monitor_vs_OPM.py which reads in the informations stored in results/Monitor_PMT_volts_summary.txt and results/(new_)pulse_power_pd_summary.txt. This code first matches the entries between the two files based on the laser pusle width and then plots the monitor PMT mean integrated value with the sigma as an error as a function of the OPM value and fits a straight line. It also plots the monitor PMT value as a function of the PD current and fits that with a straigth line. Both plots are saved in the plots/summary folder.
