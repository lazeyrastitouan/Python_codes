from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import sys

### Plot SB histograms combining 3 los from MockOb outputs (needs to have been transformed to fits by MockOb2fits.py or plotMockOb.py) at various resolution ### 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

parser = argparse.ArgumentParser(description="Plot SB histograms combining 3 los from MockOb outputs (needs to have been transformed to fits by MockOb2fits.py or plotMockOb.py) at various resolution and ratio between resolution. Argument: < --map (optional) mean or std if want to do the histogram of mean or stddev in squares of pixels (str) >, < --npix if map=mean or std, the number of pixels per square side (int) >")
parser.add_argument('--map',nargs=1, type=str,help='optional) mean or std if want to do the histogram of mean or std in squares of pixels (str)')
parser.add_argument('--npix',nargs=1, type=int,help='if map=mean or std, the number of pixels per square side (int)')
args    = parser.parse_args()
if args.map is None:
    suffix = ''
else:
    mapsty  = args.map[0]
    if args.npix is None:
        sys.exit("Please provide pixels per square side")
    nxn     = args.npix[0]
    suffix  = '_'+mapsty+'_'+str(nxn)+'x'+str(nxn)

red     = 4.
nbin    = 22
edges   = linspace(-21.,-16.6,nbin)
midbins = array([(edges[x]+edges[x+1])/2. for x in range(nbin-1)])
root    = "/data/titouan/Eagle_cosmo_gizmo/200_512/"
halo    = 'Halo11/sftresh_cool_low_fb_boost/snapdir_000/4cMpc/'
foutr   = "/data/titouan/Eagle_cosmo_gizmo/200_512/SB_Histogram_ratio_res_CDM_WDM1p5_hydro_fire2_lvl7_Halo11_z"+str(int(red))+suffix+".pdf"
fout    = "/data/titouan/Eagle_cosmo_gizmo/200_512/SB_Histogram_CDM_WDM1p5_hydro_fire2_lvl7_Halo11_z"+str(int(red))+suffix+".pdf"
los     = ['X','Y','Z']
zoom    = ['1p5','2','2p5']
DM      = ['CDM','WDM_1p5keV']
colors  = ['r','g'] 
lsty    = [':','--','-']

#prepare figure
plt.figure(1)
plt.xlim(min(midbins),max(midbins))
plt.ylim(0.5,4.0)
plt.xlabel(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$ (1+z)$^{-4}$]',fontsize=18)
plt.ylabel(r'pixel count',fontsize=18)
#for legends
plt.plot([1e-35,1e-34],[1,1],'r-',label=r'CDM')
plt.plot([1e-35,1e-34],[1,1],'g-',label=r'WDM 1.5 keV')
plt.plot([1e-35,1e-34],[1,1],'k:',label=r'$m_b = 6.7 \cdot 10^5 M_\odot$')
plt.plot([1e-35,1e-34],[1,1],'k--',label=r'$m_b = 1.0 \cdot 10^5 M_\odot$')
plt.plot([1e-35,1e-34],[1,1],'k-',label=r'$m_b = 1.6 \cdot 10^4 M_\odot$')
plt.legend(loc=0,prop={'size':14})

meanall, stdall = [[] for x in range(len(DM))],[[] for x in range(len(DM))]
for z in zoom:
    iz = zoom.index(z)
    for m in DM:
        im = DM.index(m)
        #will be removed when sim ready
        if m == 'WDM_1p5keV' and z == '2p5': 
            continue
        else:
            histo = []
            #get 3 los SB
            for l in los:
                path = root+m+'/zoom'+z+'_hydro_fire2_lvl7/'+halo+'mockCube'+l+'_s1_n0_SB_intrinsic'+suffix+'.fits'
                dhu=fits.open(path)
                SB = dhu[0].data
                SB  = array(SB)/(1.+red)**4.
                SBf = SB.flatten()
                #get rid of zeros in case noise was set to 0
                SBf[where(SBf==0.0)]=1e-35

                #get histogram
                hist   = histogram(log10(SBf),bins=edges) #array of counts (hist) and bin edges (len(hist)+1)
                histo.append(hist[0])

                del SB, SBf, hist

            #get mean and variance of 3 los
            histo = array(histo)
            meanh = mean(histo,axis=0)
            stdh  = std(histo,axis=0)
            meanall[im].append(meanh)
            stdall[im].append(stdh)
            
            #plot
            plt.plot(midbins,log10(meanh),ls=lsty[iz],color=colors[im])
            plt.fill_between(midbins,log10(meanh-stdh),log10(meanh+stdh),color=colors[im],alpha=0.2)

            del histo, meanh, stdh

plt.savefig(fout, bbox_inches='tight')
print("Saved figure to "+fout)

#plot ratio to low resolution
meanall = array(meanall)
stdall  = array(stdall)
plt.figure(2)
plt.xlim(min(midbins),max(midbins))
plt.ylim(0.0,2.5)
plt.xlabel(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$ (1+z)$^{-4}$]',fontsize=18)
plt.ylabel(r'ratio pixel count',fontsize=18)
#for legends
plt.plot([1e-35,1e-34],[1,1],'r-',label=r'CDM')
plt.plot([1e-35,1e-34],[1,1],'g-',label=r'WDM 1.5 keV')
plt.plot([1e-35,1e-34],[1,1],'k:',label=r'$m_b = 6.7 \cdot 10^5 M_\odot$')
plt.plot([1e-35,1e-34],[1,1],'k--',label=r'$m_b = 1.0 \cdot 10^5 M_\odot$')
#plt.plot([1e-35,1e-34],[1,1],'k-',label=r'$m_b = 1.6 \cdot 10^4 M_\odot$')
plt.legend(loc=0,prop={'size':14})
for im in range(len(DM)):
    for iz in range(len(zoom)):
        if iz == 2:
            continue
        if iz == 0:
            stdratio = 0.
        else:
            stdratio = abs(meanall[im,iz]/meanall[im,0])*sqrt((stdall[im,iz]/meanall[im,iz])**2.+(stdall[im,0]/meanall[im,0])**2.)
        plt.plot(midbins,meanall[im,iz]/meanall[im,0],ls=lsty[iz],color=colors[im])
        plt.fill_between(midbins,meanall[im,iz]/meanall[im,0]-stdratio,meanall[im,iz]/meanall[im,0]+stdratio,color=colors[im],alpha=0.2)

plt.savefig(foutr, bbox_inches='tight')
print("Saved figure to "+foutr)
#plt.show()
