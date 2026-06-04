from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse
import matplotlib.ticker as mticker

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

# ---
root     = '/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo'
MAGGz4   = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_4.1.npz") 
MAGGz3p5 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_3.5.npz")
noise    = '5e-20'

#prepare plot
z      = array([3.0,3.5,4.0])
zmock  = linspace(2.9,4.2)
#eagle  = (1.+zmock)**7.3/(1.+4.)**7.3
theo   = (1.+zmock)**5./(1.+4)**5.
fout   = '/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/zevol_SBprofile_dalya_zoom1p5_HELLO_MAGG_Mlow.pdf'
fig    = plt.figure()
ax     = fig.add_subplot(111)
plt.yscale('log')
plt.plot([2.9,4.2],[1,1],'k-')
#plt.plot(zmock,eagle,'k:',marker='',label=r'Eagle')
plt.plot(zmock,theo,'k--',marker='',label=r'$(1+z)^5$')
#ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
colors = ['cyan','blue']
M      = ['12.3','12.6'] 

#--- HELLO
HLOz3p5 = load('/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_z_3.596_s1_n'+noise+'_Mmin_12.5_Mhigh_13.0.npz')
HLOz4   = load('/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_z_4.022_s1_n'+noise+'_Mmin_12.5_Mhigh_13.0.npz')
filterH = HLOz3p5['binmean']<400. #ckpc, would correspond to 100pkpc at z=3
rHLOL, p16HL, p84HL = [],[],[]
for l in ['X','Y','Z']:
	rHLOtab = array([HLOz3p5['sbprofile'+l][filterH]/HLOz4['sbprofile'+l][filterH]*4.596**4./5.022**4,HLOz4['sbprofile'+l][filterH]/HLOz4['sbprofile'+l][filterH]])
	rHLOL.append(nanmedian(rHLOtab,axis=1))
	p16HL.append(nanpercentile(rHLOtab,16,axis=1))
	p84HL.append(nanpercentile(rHLOtab,84,axis=1))
rHLO = nanmedian(rHLOL,axis=0)
p16H = nanmedian(p16HL,axis=0)
p84H = nanmedian(p84HL,axis=0)
#plot
plt.plot([3.596,4.022],rHLO,'go',label=r'HELLO, $\log M = 12.9$')
plt.errorbar([3.596,4.022],rHLO,xerr=0.0,yerr=[rHLO-p16H,p84H-rHLO],fmt='go')

#--- MAGG
MAGGz4   = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_4.1.npz")
MAGGz3p5 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_3.5.npz")
rMAGGtab = array([MAGGz3p5['sbprofile']/MAGGz4['sbprofile']*4.5**4./5.1**4.,MAGGz4['sbprofile']/MAGGz4['sbprofile']])
rMAGG    = nanmedian(rMAGGtab,axis=1)
p16MAGG  = nanpercentile(rMAGGtab,16,axis=1)
p84MAGG  = nanpercentile(rMAGGtab,84,axis=1)

#plot
plt.plot([3.53,4.1],rMAGG,'ko',label=r'MAGG')
plt.errorbar([3.53,4.1],rMAGG,xerr=0.0,yerr=[rMAGG-p16MAGG,p84MAGG-rMAGG],fmt='ko')

#--DaLya
i=0
for mlow in [12.0,12.5]:
	#LOAD
	if mlow == 12.5:
		fCDMz3   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_3.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.6)+'.npz'
		fWDM3p5z3   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_3p5keV_observed_cubex_neb_z_3.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.6)+'.npz'
		fWDM1p5z3   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_1p5keV_observed_cubex_neb_z_3.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.6)+'.npz'
	else:
		fCDMz3   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_3.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
		fWDM3p5z3   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_3p5keV_observed_cubex_neb_z_3.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
		fWDM1p5z3   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_1p5keV_observed_cubex_neb_z_3.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'

	fCDMz3p5 = 'SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_3.5_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
	fCDMz4   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_4.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
	CDMz3    = load(root+'/'+fCDMz3)
	CDMz3p5  = load(root+'/'+fCDMz3p5)
	CDMz4    = load(root+'/'+fCDMz4)

	fWDM3p5z3p5 = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_3p5keV_observed_cubex_neb_z_3.5_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
	fWDM3p5z4   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_3p5keV_observed_cubex_neb_z_4.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
	WDM3p5z3    = load(root+'/'+fWDM3p5z3)
	WDM3p5z3p5  = load(root+'/'+fWDM3p5z3p5)
	WDM3p5z4    = load(root+'/'+fWDM3p5z4)

	fWDM1p5z3p5 = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_1p5keV_observed_cubex_neb_z_3.5_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
	fWDM1p5z4   = 'SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_1p5keV_observed_cubex_neb_z_4.0_s1_n'+noise+'_Mmin_'+str(mlow)+'_Mhigh_'+str(mlow+0.5)+'.npz'
	WDM1p5z3    = load(root+'/'+fWDM1p5z3)
	WDM1p5z3p5  = load(root+'/'+fWDM1p5z3p5)
	WDM1p5z4    = load(root+'/'+fWDM1p5z4)


	#GET RATIOS
	filter  = CDMz3['binmean']<400. #to consider only up to r=100pkpc at z=3
	if mlow == 12.5:
		CDMz4 = nanmedian([CDMz4['sbprofileX'][filter],CDMz4['sbprofileY'][filter],CDMz4['sbprofileZ'][filter]],axis=0)
		rCDMtab = array([CDMz3['sbprofile'][filter]/CDMz4*256./625.,CDMz3p5['sbprofile'][filter]/CDMz4*4.5**4./625.,CDMz4/CDMz4])
	else:
		rCDMtab = array([CDMz3['sbprofile'][filter]/CDMz4['sbprofile'][filter]*256./625.,CDMz3p5['sbprofile'][filter]/CDMz4['sbprofile'][filter]*4.5**4./625.,CDMz4['sbprofile'][filter]/CDMz4['sbprofile'][filter]])
	rCDM    = nanmedian(rCDMtab,axis=1)
	p16CDM  = nanpercentile(rCDMtab,16,axis=1)
	p84CDM  = nanpercentile(rCDMtab,84,axis=1)

	#rWDM3p5tab = array([WDM3p5z3['sbprofile'][filter]/WDM3p5z4['sbprofile'][filter]*256./625.,WDM3p5z3p5['sbprofile'][filter]/WDM3p5z4['sbprofile'][filter]*4.5**4./625.,WDM3p5z4['sbprofile'][filter]/WDM3p5z4['sbprofile'][filter]])
	#rWDM3p5    = nanmedian(rWDM3p5tab,axis=1)
	#p16WDM3p5  = nanpercentile(rWDM3p5tab,16,axis=1)
	#p84WDM3p5  = nanpercentile(rWDM3p5tab,84,axis=1)

	#rWDM1p5tab = array([WDM1p5z3['sbprofile'][filter]/WDM1p5z4['sbprofile'][filter]*256./625.,WDM1p5z3p5['sbprofile'][filter]/WDM1p5z4['sbprofile'][filter]*4.5**4./625.,WDM1p5z4['sbprofile'][filter]/WDM1p5z4['sbprofile'][filter]])
	#rWDM1p5    = nanmedian(rWDM1p5tab,axis=1)
	#p16WDM1p5  = nanpercentile(rWDM1p5tab,16,axis=1)
	#p84WDM1p5  = nanpercentile(rWDM1p5tab,84,axis=1)


	#PLOT
	plt.plot(z,rCDM,marker='o',color=colors[i],ls='',label=r'DaLya, $\log M =$ '+M[i])
	plt.errorbar(z,rCDM,xerr=0.0,yerr=[rCDM-p16CDM,p84CDM-rCDM],marker='o',color=colors[i],ls='')
	#plt.plot(z+0.015,rWDM3p5,'go',label=r'WDM 3.5 keV')
	#plt.errorbar(z+0.015,rWDM3p5,xerr=0.0,yerr=[rWDM3p5-p16WDM3p5,p84WDM3p5-rWDM3p5],fmt='go')
	#plt.plot(z-0.015,rWDM1p5,'ro',label=r'WDM 1.5 keV')
	#plt.errorbar(z-0.015,rWDM1p5,xerr=0.0,yerr=[rWDM1p5-p16WDM1p5,p84WDM1p5-rWDM1p5],fmt='ro')

	i+=1

#finalize plot
plt.xlim(2.9,4.2)
plt.ylim([0.25,2.0])
plt.xlabel(r'z',fontsize=14)
plt.ylabel(r'$[SB(z)\cdot(1+z)^4]/[SB(z=4)\cdot(1+4)^4]$', fontsize=14)
plt.legend(loc=0,prop={'size':12},frameon=False)
plt.savefig(fout, bbox_inches='tight')
plt.show()
