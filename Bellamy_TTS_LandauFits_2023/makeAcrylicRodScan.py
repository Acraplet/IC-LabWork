##This is a little python code that extracts the information in the .txt file output of the Bellamy analysis and plots the results to make scans. This is a 1D scan for now

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


list_names = ['runNumber', 'channelTrigger', 'laserPower', 'PMT1_number', 'PMT1_theta', 'PMT2_number', 'PMT2_theta', 'PMT3_number', 'PMT3_theta', 'PMT4_name', 'PMT4_theta', 'LBphi', 'PMT1_gain', 'PMT2_gain', 'PMT3_gain', 'PMT4_gain', 'PMT1_gain_err', 'PMT2_gain_err', 'PMT3_gain_err', 'PMT4_gain_err', 'PMT1_Q1', 'PMT2_Q1', 'PMT3_Q1', 'PMT4_Q1', 'PMT1_1peBel', 'PMT2_1peBel', 'PMT3_1peBel', 'PMT4_1peBel', 'PMT1_2peBel', 'PMT2_2peBel', 'PMT3_2peBel', 'PMT4_2peBel', 'PMT1_chi2perPt', 'PMT2_chi2perPt', 'PMT3_chi2perPt', 'PMT4_chi2perPt', 'totalDataDuration', 'PMT1_nHits_TTS', 'PMT2_nHits_TTS', 'PMT3_nHits_TTS', 'PMT4_nHits_TTS', 'PMT1_nHits_TTS_err', 'PMT2_nHits_TTS_err', 'PMT3_nHits_TTS_err', 'PMT4_nHits_TTS_err', 'PMT1_nHits_DR', 'PMT2_nHits_DR', 'PMT3_nHits_DR', 'PMT4_nHits_DR', 'PMT1_nHits_DR_err', 'PMT2_nHits_DR_err', 'PMT3_nHits_DR_err', 'PMT4_nHits_DR_err', 'PMT1_FWHM', 'PMT2_FWHM', 'PMT3_FWHM', 'PMT4_FWHM', 'PMT1_FWHM_err', 'PMT2_FWHM_err', 'PMT3_FWHM_err', 'PMT4_FWHM_err']

PMT_names = ['Trigger', 'Square', 'Circle', 'Losange', 'Triangle', 'Monitor']



df = pd.read_csv('acrylicRodScan.txt', sep=' ', header = None, names = list_names)

df = df.sort_values(by="LBphi")

df = df[df["laserPower"]==130]
# df = df[df["LBphi"]==183]
df = df[df["runNumber"]<=90]

marker = ['x', 'x']
colors = ['blue', 'orange', 'lightgray']

print(df, print(df["PMT2_nHits_TTS"]/(df["PMT3_nHits_DR"]+df["PMT3_nHits_TTS"])))

for PMT in [1, 2, 3]:
        #print()
        for powerID in range(len(df["laserPower"].unique())):
            power = df["laserPower"].unique()[powerID]
            df_buf = df[df["laserPower"]==power]


            totalMonitorHits = df["PMT3_nHits_DR"]+df["PMT3_nHits_TTS"]


            onePE = df_buf['PMT%i_nHits_TTS'%(PMT)]/(totalMonitorHits)
            onePE_err =  np.sqrt(onePE*totalMonitorHits)/totalMonitorHits #/(totalMonitorHits) df_buf['PMT%i_1peBel_err'%(PMT)]
            plt.errorbar(df_buf['runNumber'], onePE, xerr = 0.05, yerr=onePE_err, label = 'PMT%i: %s, number of PE/number of PE in Monitor, Laser %imA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], power), fmt = '%s'%(marker[powerID]), color = colors[powerID])



        plt.title("Number of TTS fitted PE - 120mA, 183 degrees")
        plt.xlabel('Run Number')
        plt.ylabel('Number of PE per PE in Monitor PMT')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.show()


for PMT in [1, 2, 3, 4]:
    df['PMT%s_totpeBel'%PMT] = df['PMT%i_1peBel'%PMT] + df['PMT%i_1peBel'%PMT]*2

lineType=['-', '--', '-']

for value in ["DR", "TTS"]:
    for PMT in [1, 2]:
        #print()
         for power in df["laserPower"].unique():
            df_buf = df[df["laserPower"]==power]
            plt.errorbar(df_buf['runNumber'], df_buf['PMT%i_nHits_%s'%(PMT, value)]/(df_buf['totalDataDuration']*1000), xerr = 0.05, yerr=df_buf['PMT%i_nHits_%s_err'%(PMT, value)]/(df_buf['totalDataDuration']*1000), label = 'PMT%i: %s, %s laser %i mA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], value, power), fmt = "%s"%lineType[PMT-1])
        
            
    plt.title("Number of hits in %s as a function of Acrylic rod orientation"%value)
    plt.xlabel('Run Number')
    plt.ylabel('Number of events per 1000 triggers (i.e. kHz)')
    plt.grid()
    plt.legend()
    plt.show()
    

monitorPMTNumber = 3

plt.figure(figsize = (10,20))
for value in ["TTS"]:
    for PMT in [1, 2]:
        #print()
        for power in df["laserPower"].unique():
            df_buf = df[df["laserPower"]==power]
            totalMonitorHits = df_buf['PMT%i_nHits_TTS'%(monitorPMTNumber)] + df_buf['PMT%i_nHits_DR'%(monitorPMTNumber)]
            plt.errorbar(df_buf['LBphi']-183, df_buf['PMT%i_nHits_%s'%(PMT, value)]/(totalMonitorHits), xerr = 0.05, yerr=df_buf['PMT%i_nHits_%s_err'%(PMT, value)]/(totalMonitorHits), label = 'PMT%i: %idegrees\n #signal hits in %s over #hits in Monitor PMT \n Laser %imA'%(PMT, np.array(df_buf["PMT%i_theta"%PMT])[0]-90, value, power), fmt = 'x--')
        
            
    plt.title("Acrylic light guide characterisation\nNumber of signal hits in %s corrected by monitor PMT signal"%value)
    plt.xlabel('Angular orienation (degrees)')
    plt.ylabel('Number of signal events per monitor PMT hit')
    # plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.savefig("AcrylicRodScan.pdf")
    plt.show()


colors = ['blue', 'orange', 'lightgray']
marker = ['x', 'x']
for PMT in [1, 2, 3]:
        #print()
        for powerID in range(len(df["laserPower"].unique())):
            power = df["laserPower"].unique()[powerID]
            df_buf = df[df["laserPower"]==power]
            onePE = df_buf['PMT%i_1peBel'%(PMT)]/(df_buf['totalDataDuration']*1000) #/(totalMonitorHits)
            onePE_err =  np.sqrt(onePE)/(df_buf['totalDataDuration']*1000) #/(totalMonitorHits) df_buf['PMT%i_1peBel_err'%(PMT)]
            plt.errorbar(df_buf['LBphi'], onePE, xerr = 0.05, yerr=onePE_err, label = 'PMT%i: %s, number of 1PE events (bellamy), Laser %imA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], power), fmt = '%s-'%(marker[powerID]), color = colors[powerID])

            twoPE = df_buf['PMT%i_2peBel'%(PMT)]/(df_buf['totalDataDuration']*1000) #/(totalMonitorHits)
            twoPE_err =  np.sqrt(twoPE)/(df_buf['totalDataDuration']*1000) #/(totalMonitorHits) df_buf['PMT%i_1peBel_err'%(PMT)]
            plt.errorbar(df_buf['LBphi'], twoPE, xerr = 0.05, yerr=twoPE_err, label = 'PMT%i: %s, number of 2PE events (bellamy), Laser %imA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], power), fmt = '%s--'%(marker[powerID]), color = colors[powerID])

            plt.errorbar(df_buf['LBphi'], twoPE * 2 + onePE, xerr = 0.05, yerr=np.sqrt(2*twoPE_err**2 + onePE_err**2), label = 'PMT%i: %s, number of photons (1pe+2pe bellamy), Laser %imA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], power), fmt = '%s:'%(marker[powerID]), color = colors[powerID])

        plt.title("Number of Bellamy fitted PE as a function of Acrylic rod orientation")
        plt.xlabel('Orientation (degrees)')
        plt.ylabel('Number of events per 1000 triggers (i.e. kHz)')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.show()

for PMT in [1, 2, 3]:
        #print()
        for powerID in range(len(df["laserPower"].unique())):
            power = df["laserPower"].unique()[powerID]
            df_buf = df[df["laserPower"]==power]


            totalMonitorHits = df_buf['PMT%i_totpeBel'%(monitorPMTNumber)]


            onePE = df_buf['PMT%i_totpeBel'%(PMT)]/(totalMonitorHits)
            onePE_err =  np.sqrt(onePE*totalMonitorHits)/totalMonitorHits #/(totalMonitorHits) df_buf['PMT%i_1peBel_err'%(PMT)]
            plt.errorbar(df_buf['LBphi'], onePE, xerr = 0.05, yerr=onePE_err, label = 'PMT%i: %s, number of PE/number of PE in Monitor, Laser %imA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], power), fmt = '%s-'%(marker[powerID]), color = colors[powerID])



        plt.title("Number of Bellamy fitted PE as a function of Acrylic rod orientation")
        plt.xlabel('Orientation (degrees)')
        plt.ylabel('Number of PE per PE in Monitor PMT')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.show()



