##This is a little python code that extracts the information in the .txt file output of the Bellamy analysis and plots the results to make scans. This is a 1D scan for now

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT


list_names = ['runNumber', 'channelTrigger', 'laserPower', 'PMT1_number', 'PMT1_theta', 'PMT2_number', 'PMT2_theta', 'PMT3_number', 'PMT3_theta', 'PMT4_name', 'PMT4_theta', 'LBphi', 'PMT1_gain', 'PMT2_gain', 'PMT3_gain', 'PMT4_gain', 'PMT1_gain_err', 'PMT2_gain_err', 'PMT3_gain_err', 'PMT4_gain_err', 'PMT1_Q1', 'PMT2_Q1', 'PMT3_Q1', 'PMT4_Q1', 'PMT1_1peBel', 'PMT2_1peBel', 'PMT3_1peBel', 'PMT4_1peBel', 'PMT1_2peBel', 'PMT2_2peBel', 'PMT3_2peBel', 'PMT4_2peBel', 'PMT1_chi2perPt', 'PMT2_chi2perPt', 'PMT3_chi2perPt', 'PMT4_chi2perPt', 'totalDataDuration', 'PMT1_nHits_TTS', 'PMT2_nHits_TTS', 'PMT3_nHits_TTS', 'PMT4_nHits_TTS', 'PMT1_nHits_TTS_err', 'PMT2_nHits_TTS_err', 'PMT3_nHits_TTS_err', 'PMT4_nHits_TTS_err', 'PMT1_nHits_DR', 'PMT2_nHits_DR', 'PMT3_nHits_DR', 'PMT4_nHits_DR', 'PMT1_nHits_DR_err', 'PMT2_nHits_DR_err', 'PMT3_nHits_DR_err', 'PMT4_nHits_DR_err', 'PMT1_FWHM', 'PMT2_FWHM', 'PMT3_FWHM', 'PMT4_FWHM', 'PMT1_FWHM_err', 'PMT2_FWHM_err', 'PMT3_FWHM_err', 'PMT4_FWHM_err']

PMT_names = ['Trigger', 'Square', 'Circle', 'Losange', 'Triangle', 'Monitor']



df = pd.read_csv('LaserBallPrototype1Scan.txt', sep=' ', header = None, names = list_names)

df = df.sort_values(by="runNumber")

Monitor = 3

df = df[df['laserPower']>= 1] #get rid of DR
df = df[df['runNumber']>=109]

#The error on the angular resolution 0.05 for the phi, 2.5 for the theta
x_width = 2.5
x_width2 = 0.05

#Decide if we want the monitor PMT to also have the Dr in
useMonitorDR = True
phiCount = 0

colors = ['black', 'red', 'green', 'blue', 'orange', 'gray']
#['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
 #'#7f7f7f', '#bcbd22', '#17becf']
#
plt.figure(1, figsize=(8, 6))

target = "nHits_TTS"
target_error = "nHits_TTS_err"

title = "Scan of Laser Ball, photon escape time spread"
ylabel = "TTS peak FWHM(ns)"
xlabel='Theta (degrees)'
#"Scan of Laser Ball, stability of the dark rate"

for i, phi in enumerate(np.sort(df['LBphi'].unique())):

    df_buf_phi = df[df['LBphi'] == phi]
    for power in df_buf_phi['laserPower'].unique():
        count = 0
        df_buf = df_buf_phi[df_buf_phi['laserPower'] == power]
        #calculate the monitor PMT signal, with or without DR depending
        if useMonitorDR == True:
            nphotonsMonitor= df_buf["PMT%i_nHits_TTS"%Monitor] + df_buf["PMT%i_nHits_DR"%Monitor]
            monitor_err = np.sqrt(df_buf["PMT%i_nHits_TTS_err"%Monitor]**2 +  df_buf["PMT%i_nHits_DR_err"%Monitor]**2 )
        else:
            nphotonsMonitor=df_buf["PMT%i_nHits_TTS"%Monitor]
            monitor_err = df_buf["PMT%i_nHits_TTS_err"%Monitor]

        for PMT in [1,2]:
            nphotonsPMT = df_buf["PMT%i_%s"%(PMT, target)]
            PMT_err = df_buf["PMT%i_%s"%(PMT, target_error)]
            ratio = nphotonsPMT/nphotonsMonitor
            ratio_error = np.sqrt((PMT_err/nphotonsMonitor)**2 + (ratio * monitor_err / nphotonsMonitor)**2 )
            x_array = df_buf['PMT%i_theta'%PMT]
            #df_buf['PMT%i_number'%PMT]#df_buf['runNumber']#df_buf['PMT%i_theta'%PMT]

            for j in range(len(ratio)):
                ratio = np.array(ratio)
                ratio_error = np.array(ratio_error)

                print("The fractional statistical error for 0.25s of data at %imA at phi = %i and theta = %i is %.4f"%(power, phi, np.array(x_array)[j],ratio_error[j]/ratio[j]))

                plt.figure(1)
                plt.fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[ratio[j] - ratio_error[j], ratio[j] - ratio_error[j]], [ratio[j] + ratio_error[j],ratio[j] + ratio_error[j]], alpha=0.05,color = '%s'%colors[phiCount])

            if count == 0:
                plt.figure(1)
                plt.errorbar(x_array, ratio, xerr = x_width, yerr = ratio_error, label = 'Phi = %i deg - %imA'%(phi, power), fmt = 'x', color = '%s'%colors[phiCount])
                # plt.figure(2)
                # plt.errorbar(x_array, ratio, xerr = x_width2, yerr = ratio_error, label = 'Phi = %i deg'%phi, fmt = 'x', color = '%s'%colors[phiCount])
                count += 1
            else:
                plt.figure(1)
                plt.errorbar(x_array, ratio, xerr = x_width, yerr = ratio_error,  fmt = 'x', color = '%s'%colors[phiCount])
                count += 1
        phiCount += 1

#labels
plt.figure(1)
plt.xlabel('%s'%xlabel, fontsize = 15)
if useMonitorDR:
    # plt.ylabel('(Number of signal PE)/(Total number of hits in Monitor PMT)', fontsize = 15)
    plt.ylabel('%s'%ylabel, fontsize = 15)
    # plt.title("Scan of Laser Ball, #signal PE corrected by #monitor PE", fontsize = 15, weight = 'bold')
    plt.title("%s"%(title), fontsize = 15, weight = 'bold')
else:
    plt.ylabel('Number of signal PE/Number of signal in Monitor', fontsize = 15)
    plt.title("Scan of Laser Ball, #signal PE corrected by #signal monitor PE", fontsize = 15, weight = 'bold')
plt.legend()

plt.grid()
plt.show()

# #Make the final plots
# for PMT in [1, 2, 3, 4]:
#     for i, phi in enumerate(np.sort(df['LBphi'].unique())):
#         df_buf_phi = df[df['LBphi'] == phi]
#         plt.plot()

target = "nHits_TTS"
target_error = "nHits_TTS_err"
phiCount = 0

fig, axs = plt.subplots(2, 1, figsize=(10, 15), gridspec_kw={'height_ratios': [3, 1]})

for theta in [0, 30, 90]:
    for PMT in [1, 2, 3, 4]:

        # plt.subplot(2, 1, 1)
        # plt.gca().set_aspect('equal', adjustable='box')
        df_buf = df[df['PMT%i_theta'%PMT]==theta]
        x_array = df_buf['LBphi']

        nphotonsPMT = df_buf["PMT%i_%s"%(PMT, target)]
        PMT_err = df_buf["PMT%i_%s"%(PMT, target_error)]
        nphotonsMonitor=df_buf["PMT%i_nHits_TTS"%Monitor] + df_buf["PMT%i_nHits_DR"%Monitor]
        monitor_err = df_buf["PMT%i_nHits_TTS_err"%Monitor]

        ratio = nphotonsPMT/nphotonsMonitor
        ratio_error = np.sqrt((PMT_err/nphotonsMonitor)**2 + (ratio * monitor_err / nphotonsMonitor)**2 )

        ratio = np.array(ratio)
        ratio_error = np.array(ratio_error)



        for j in range(len(ratio)):

            axs[0].fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[ratio[j] - ratio_error[j], ratio[j] - ratio_error[j]], [ratio[j] + ratio_error[j],ratio[j] + ratio_error[j]], alpha=0.05, color = '%s'%colors[phiCount])

        if len(ratio)!=0:

            axs[0].errorbar(df_buf['LBphi'], ratio, yerr=ratio_error, fmt = 'x', color = "%s"%colors[phiCount])
            axs[0].plot([x_array.min()*0.95, x_array.max()*1.05], [ratio.mean(), ratio.mean()], color = '%s'%colors[phiCount], linestyle='--', label = 'Theta = %i degrees Mean: %.2e +/- %.2e (%.2f percent deviation)'%(theta, ratio.mean(),ratio.std(), ratio.std()/ratio.mean() * 100))


        if len(ratio)!=0:

            axs[1].errorbar(df_buf['LBphi'], (ratio-ratio.mean())/ratio.mean(), yerr=ratio_error/ratio.mean(), fmt = 'x', color = "%s"%colors[phiCount])


            axs[1].plot([x_array.min()*0.95, x_array.max()*1.05], [0, 0], color = 'darkgray', linestyle='--')

        for j in range(len(ratio)):

            axs[1].fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[(ratio[j]-ratio.mean())/ratio.mean() - ratio_error[j]/ratio.mean(), (ratio[j]-ratio.mean())/ratio.mean() - ratio_error[j]/ratio.mean()], [(ratio[j]-ratio.mean())/ratio.mean() + ratio_error[j]/ratio.mean(), (ratio[j]-ratio.mean())/ratio.mean() + ratio_error[j]/ratio.mean()], alpha=0.05, color = '%s'%colors[phiCount])




        # plt.fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[ratio[j] - ratio_error[j], ratio[j] - ratio_error[j]], [ratio[j] + ratio_error[j],ratio[j] + ratio_error[j]], alpha=0.05,color = '%s'%colors[phiCount])
    phiCount += 1

axs[0].legend()
axs[1].grid()
axs[0].grid()
axs[0].set_ylim(0.034, 0.044)
axs[1].set_xlabel('Phi (degrees)', fontsize = 14)
axs[0].set_xlabel('Phi (degrees)', fontsize = 14)
axs[0].set_title('WCTE LB prototype 1 - Amount of light (laser fluctuations corrected)', weight = 'bold', fontsize = 14)
axs[0].set_ylabel('Number of signal hits corrected by monitor PMT', fontsize = 14)
axs[1].set_ylabel('data-mean/data', fontsize = 14)
# plt.savefig("prototype1_characterisation.pdf")
plt.show()

target = "FWHM"
target_error = "FWHM_err"
phiCount = 0

stdToFWHM = 1#2.355
PMT_TTS = [0.4075, 0.4192, 0.4136, 0.4153]

fig, axs = plt.subplots(2, 1, figsize=(10, 15), gridspec_kw={'height_ratios': [3, 1]})

for theta in [0, 30, 90, 120]:
    for PMT in [1, 2, 3, 4]:

        # plt.subplot(2, 1, 1)
        # plt.gca().set_aspect('equal', adjustable='box')
        df_buf = df[df['PMT%i_theta'%PMT]==theta]
        x_array = df_buf['LBphi']

        nphotonsPMT = df_buf["PMT%i_%s"%(PMT, target)]
        print(nphotonsPMT)
        PMT_err = df_buf["PMT%i_%s"%(PMT, target_error)]
        nphotonsMonitor=1 ##df_buf["PMT%i_nHits_TTS"%Monitor] + df_buf["PMT%i_nHits_DR"%Monitor]
        monitor_err = 0 #df_buf["PMT%i_nHits_TTS_err"%Monitor]

        if len(nphotonsPMT)!=0:
            nphotonsPMT = nphotonsPMT - stdToFWHM * PMT_TTS[int(df_buf["PMT%i_number"%PMT].unique()[0])-1]

        ratio = nphotonsPMT/nphotonsMonitor
        ratio_error = np.sqrt((PMT_err/nphotonsMonitor)**2 + (ratio * monitor_err / nphotonsMonitor)**2 )

        ratio = np.array(ratio)
        ratio_error = np.array(ratio_error)



        for j in range(len(ratio)):


            axs[0].fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[ratio[j] - ratio_error[j], ratio[j] - ratio_error[j]], [ratio[j] + ratio_error[j],ratio[j] + ratio_error[j]], alpha=0.05, color = '%s'%colors[phiCount])

        if len(ratio)!=0:

            axs[0].errorbar(df_buf['LBphi'], ratio, yerr=ratio_error, fmt = 'x', color = "%s"%colors[phiCount])
            axs[0].plot([x_array.min()*0.95, x_array.max()*1.05], [ratio.mean(), ratio.mean()], color = '%s'%colors[phiCount], linestyle='--', label = 'PMT %i: Mean: %.2e +/- %.2e (%.2f percent deviation)'%(df_buf["PMT%i_number"%PMT].unique()[0], ratio.mean(),ratio.std(), ratio.std()/ratio.mean() * 100))


        if len(ratio)!=0:

            axs[1].errorbar(df_buf['LBphi'], (ratio-ratio.mean())/ratio.mean(), yerr=ratio_error/ratio.mean(), fmt = 'x', color = "%s"%colors[phiCount])


            axs[1].plot([x_array.min()*0.95, x_array.max()*1.05], [0, 0], color = 'darkgray', linestyle='--')

        for j in range(len(ratio)):

            axs[1].fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[(ratio[j]-ratio.mean())/ratio.mean() - ratio_error[j]/ratio.mean(), (ratio[j]-ratio.mean())/ratio.mean() - ratio_error[j]/ratio.mean()], [(ratio[j]-ratio.mean())/ratio.mean() + ratio_error[j]/ratio.mean(), (ratio[j]-ratio.mean())/ratio.mean() + ratio_error[j]/ratio.mean()], alpha=0.05, color = '%s'%colors[phiCount])




        # plt.fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[ratio[j] - ratio_error[j], ratio[j] - ratio_error[j]], [ratio[j] + ratio_error[j],ratio[j] + ratio_error[j]], alpha=0.05,color = '%s'%colors[phiCount])
    phiCount += 1

axs[0].legend()
axs[1].grid()
axs[0].grid()
# axs[0].set_ylim(0.65, 1)
axs[1].set_xlabel('Phi (degrees)', fontsize = 14)
axs[0].set_xlabel('Phi (degrees)', fontsize = 14)
axs[0].set_title('WCTE LB prototype 1 - Photon escape time spread', weight = 'bold', fontsize = 14)
axs[0].set_ylabel('Laser Ball photon escape time (ns)', fontsize = 14)
axs[1].set_ylabel('data-mean/data', fontsize = 14)
# plt.savefig("prototype1_characterisation.pdf")
plt.show()





raise end

#
#     df_buf_phi = df[df['LBphi'] == phi]
#     for power in df_buf_phi['laserPower'].unique():
#         count = 0
#         df_buf = df_buf_phi[df_buf_phi['laserPower'] == power]
#         #calculate the monitor PMT signal, with or without DR depending
#         if useMonitorDR == True:
#             nphotonsMonitor= 1#df_buf["PMT%i_nHits_TTS"%Monitor] + df_buf["PMT%i_nHits_DR"%Monitor]
#             monitor_err = np.sqrt(df_buf["PMT%i_nHits_TTS_err"%Monitor]**2 +  df_buf["PMT%i_nHits_DR_err"%Monitor]**2 )
#         else:
#             nphotonsMonitor=df_buf["PMT%i_nHits_TTS"%Monitor]
#             monitor_err = df_buf["PMT%i_nHits_TTS_err"%Monitor]
#
#         for PMT in [1,2]:
#             nphotonsPMT = df_buf["PMT%i_%s"%(PMT, target)]
#             PMT_err = df_buf["PMT%i_%s"%(PMT, target_error)]
#             ratio = nphotonsPMT#nphotonsPMT/nphotonsMonitor
#             ratio_error = PMT_err#np.sqrt((PMT_err/nphotonsMonitor)**2 + (ratio * monitor_err / nphotonsMonitor)**2 )
#             x_array = df_buf['PMT%i_theta'%PMT]
#             #df_buf['PMT%i_number'%PMT]#df_buf['runNumber']#df_buf['PMT%i_theta'%PMT]
#
#             for j in range(len(ratio)):
#                 ratio = np.array(ratio)
#                 ratio_error = np.array(ratio_error)
#
#                 print("The fractional statistical error for 0.25s of data at %imA at phi = %i and theta = %i is %.4f"%(power, phi, np.array(x_array)[j],ratio_error[j]/ratio[j]))
#
#                 plt.figure(1)
#                 plt.fill_between([np.array(x_array)[j]-x_width, np.array(x_array)[j]+x_width],[ratio[j] - ratio_error[j], ratio[j] - ratio_error[j]], [ratio[j] + ratio_error[j],ratio[j] + ratio_error[j]], alpha=0.05,color = '%s'%colors[phiCount])
#
#             if count == 0:
#                 plt.figure(1)
#                 plt.errorbar(x_array, ratio, xerr = x_width, yerr = ratio_error, label = 'Phi = %i deg - %imA'%(phi, power), fmt = 'x', color = '%s'%colors[phiCount])
#                 # plt.figure(2)
#                 # plt.errorbar(x_array, ratio, xerr = x_width2, yerr = ratio_error, label = 'Phi = %i deg'%phi, fmt = 'x', color = '%s'%colors[phiCount])
#                 count += 1
#             else:
#                 plt.figure(1)
#                 plt.errorbar(x_array, ratio, xerr = x_width, yerr = ratio_error,  fmt = 'x', color = '%s'%colors[phiCount])
#                 count += 1
#         phiCount += 1
#
# #labels
# plt.figure(1)
# plt.xlabel('%s'%xlabel, fontsize = 15)
# if useMonitorDR:
#     # plt.ylabel('(Number of signal PE)/(Total number of hits in Monitor PMT)', fontsize = 15)
#     plt.ylabel('%s'%ylabel, fontsize = 15)
#     # plt.title("Scan of Laser Ball, #signal PE corrected by #monitor PE", fontsize = 15, weight = 'bold')
#     plt.title("%s"%(title), fontsize = 15, weight = 'bold')
# else:
#     plt.ylabel('Number of signal PE/Number of signal in Monitor', fontsize = 15)
#     plt.title("Scan of Laser Ball, #signal PE corrected by #signal monitor PE", fontsize = 15, weight = 'bold')
# plt.legend()
#
# plt.grid()
# plt.show()







raise end

# minPhi = min(df['LBphi'])
# maxPhi = max(df['LBphi'])
# minTheta = min(min(df['PMT2_theta']), min(df['PMT1_theta']))
# maxTheta = max(max(df['PMT2_theta']), max(df['PMT1_theta']))
# X, Y, Z = [], [], []


# num_bins = 4
# phi_bins = np.sort(df["LBphi"].unique())
#
# allTheta = []
# for PMT in [1,2]:
#     for i in np.array(df['PMT%i_theta'%PMT]):
#         allTheta.append(i)

#
# theta_bins = np.sort(pd.Series(allTheta).unique()) #np.linspace(minTheta, maxTheta, num_bins + 1)
# print(minPhi, maxPhi, minTheta, maxTheta)


'''
        print(df_buf["PMT%i_theta"%PMT])
        for theta in df_buf["PMT%i_theta"%PMT]:
            df_buf2 = df_buf[df_buf["PMT%i_theta"%PMT]==theta]
            nHitsCorrected =df_buf2["PMT%s_nHits_TTS"%PMT]/(df_buf2["PMT%s_nHits_TTS"%Monitor])# + df_buf2["PMT%s_nHits_TTS"%Monitor])
            # df_buf2["PMT%s_nHits_TTS"%PMT]#/
            print(nHitsCorrected)
            X.append(float(phi))
            Y.append(float(df_buf2["PMT%i_theta"%PMT]))
            Z.append(float(nHitsCorrected))


plt.grid(True)
multiplicator = 1000
Z = np.array(Z) * multiplicator
plt.scatter(X, Y, c= Z, cmap='viridis', marker = 's', s = 1000)
plt.colorbar(label = 'Number of PMT signal hits (x %i) per signal monitor hits'%multiplicator)
plt.xticks(phi_bins)
plt.yticks(theta_bins)
plt.xlabel("phi (degrees)")
plt.ylabel("theta (degrees)")
for i, z in enumerate(Z):
    if z>np.mean(Z):
        plt.text(X[i], Y[i], f'{z:.2f}', ha='center', va='center', color='black')
    else:
        plt.text(X[i], Y[i], f'{z:.2f}', ha='center', va='center', color='lightgray')
# plt.imshow(Z, extent=[minPhi, maxPhi, minTheta, maxTheta], cmap='viridis')
plt.show()

# for PMT in [1,2]:
#     for i in range(num_bins):
#         for j in range(num_bins):
#             phi_lower, phi_upper = phi_bins[i], phi_bins[i + 1]
#             theta_lower, theta_upper = theta_bins[j], theta_bins[j + 1]
#
#             # Filter dataframe based on bin ranges
#             df_buf = df[(df['LBphi'] >= phi_lower) & (df['LBphi'] < phi_upper) &
#                         (df['PMT%i_theta'%PMT] >= theta_lower) & (df['PMT%i_theta'%PMT] < theta_upper)]
#
#             # Calculate Z value (mean nHitsCorrected in the bin)
#             if len(df_buf) > 0:
#                 Z_value = df_buf["PMT1_nHits_TTS"].mean()
#             else:
#                 Z_value = np.nan
#
#             # Append values for plotting
#             X.append((phi_lower + phi_upper) / 2)
#             Y.append((theta_lower + theta_upper) / 2)
#             Z.append(Z_value)
#



raise End
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


for value in ["TTS"]:
    for PMT in [2]:
        #print()
        for power in df["laserPower"].unique():
            df_buf = df[df["laserPower"]==power]
            totalMonitorHits = df_buf['PMT%i_nHits_TTS'%(monitorPMTNumber)] + df_buf['PMT%i_nHits_DR'%(monitorPMTNumber)]
            plt.errorbar(df_buf['LBphi'], df_buf['PMT%i_nHits_%s'%(PMT, value)]/(totalMonitorHits), xerr = 0.05, yerr=df_buf['PMT%i_nHits_%s_err'%(PMT, value)]/(totalMonitorHits), label = 'PMT%i: %s/monitorPMT, %s, Laser %imA'%(PMT, PMT_names[df['PMT%i_number'%PMT].unique()[0]], value, power), fmt = 'x--')
        
            
    plt.title("Number of hits in %s as a function of Acrylic rod orientation"%value)
    plt.xlabel('Orientation (degrees)')
    plt.ylabel('Number of signal events per monitorPMT hit')
    # plt.yscale('log')
    plt.grid()
    plt.legend()
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

'''

