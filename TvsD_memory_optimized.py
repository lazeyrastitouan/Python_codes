#####################################################################################################################
#                                                                                                                   #
#       This codes plots the temperature density diagram for                                                        #
#           a) all particles                                                                                        #
#           b) particles within Rvirs of a halo catalog                                                             #
#           c) particles outside all Rvirs of a halo catalog                                                        #
#                                                                                                                   #
#       This codes does not verify if halos are main or subhalos. This distinction should be done beforhand         #
#                                                                                                                   #
#       Format of sims: gizmo/eagle hdf5 (needs keys tempeature and density in PartType0 group)                     #
#       Format of halo catalogs: AHF/P2C ASCII (at least 12 columns: Host, x, y, z, Rvir in col 2, 6, 7, 8, 12)     #  
#                                                                                                                   #
#       T. Lazeyras, 2025                                                                                           #
#                                                                                                                   #
#####################################################################################################################

from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib
#from pylab import *
#from matplotlib.ticker import ScalarFormatter
#import matplotlib.gridspec as gridspec
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import* 
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree
from scipy.stats import binned_statistic_2d

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)


def insphere(self, part, halo, radius):
    diff = part-halo
    dist = linalg.norm(diff)

    return dist <= radius

# What to do
parser = argparse.ArgumentParser(description='Plot T vs D for gizmo output high refinement. Args: < path to snapshot >, < simulation type >, < optional for inside/outside Rvir: halo catalog >, < if Hcat is provided, in or out (of Rvir) >, < optional lower and higher edge for the annulus inside Rvir (in frac of Rvir) >, < optional lower and higher edge for the annulus inside Rvir (in ckpc) >, < optional density_temperature file to plot if already saved >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('simtype',nargs=1,type=str, help='Simulation type (dalya/eagle/colibre/hello) (str, for hello rename the snapshot starting with "snap")')
parser.add_argument('--Hcat',nargs=1,type=str, help='Halo catalog (ascii), in that case plots TvsD inside/outside of halos')
parser.add_argument('--InOut',nargs=1,type=str,help='For Eagle, or in case a halo catalog is provided, if do plot for particles inside or outside (in/out)')
parser.add_argument('--Rvirfrac',nargs=2,type=float, help='Lower and higher edges of the annuli in fraction of Rvir to consider if "in" (between 0 and 1) (2 floats)')
parser.add_argument('--annulus',nargs=2,type=float,help='Lower and higher edges of the annuli in ckpc to consider if "in" (make sure higher is within CGM!) (2 floats)')
parser.add_argument('--dfile',nargs=1,type=str,help='If density and temperature have been saved to file, will load them directly from snappath/dfile (str)')
args     = parser.parse_args()
snappath = args.snappath[0]
sim      = args.simtype[0]
catbool  = False
inout    = False
annulus  = False
if args. Hcat is not None:
    catbool = True
if args.InOut is not None:
    inout   = True
    if args.InOut[0] == 'in':
        out = False
    else:
        out = True
    if args.annulus is not None:
        annulus = True

# settings
datatype  = 'HDF5' #'Gadget2'
cattype   = 'P2C' #'AHF' #'P2C'
verbose   = True
savenpy   = True
plot      = True
add_props = True
ptype     = 'PartType0'
	
print("Simulation type = "+sim)
print("Data type = "+datatype)
print("Catalog type = "+cattype)

#get density/T from simulation output
if args.dfile is None:
    h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5' and (f[0:4] == 'snap' or f[0:7]=='colibre')]
    h5files  = sorted(h5files)
    h5files  = sorted(h5files,key=len)
    numfiles = len(h5files)
    if verbose:
	    print(h5files)
    
    if sim == 'colibre':
        h = h5py.File(snappath+'/'+h5files[0],'r')['Cosmology'].attrs['h'][0]
        U_M  = float(h5py.File(snappath+'/'+h5files[0],'r')['Units'].attrs['Unit mass in cgs (U_M)'][0])
        U_L  = float(h5py.File(snappath+'/'+h5files[0],'r')['Units'].attrs['Unit length in cgs (U_L)'][0])
        U_t  = float(h5py.File(snappath+'/'+h5files[0],'r')['Units'].attrs['Unit time in cgs (U_t)'][0])
        z = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['Redshift'][0]
    else:
        h = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['HubbleParam']
        z = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['Redshift']
    a = 1./(1.+z)

    if sim == 'eagle':
        #define a region as eagle is too big
        centerH = array([24191.,77352.,61463.])*h #ckpc/h
        dx      = 20000.*h
        print("Considering particles only in "+str(centerH-dx)+" <= x <= "+str(centerH+dx))
        
    # Load particles
    dens_list, temp_list, mass_list = [], [], []
    pos_list = []
    if add_props:
        Zsol = 0.013
        Z_list, HI_list = [],[]
    if sim == 'eagle' and not catbool and inout:
        GN_list = []

    # --- HDF5 FILES LOOP --- #

    print("Loading...")
    for i in h5files:
        f   = h5py.File(snappath+"/"+i, "r")
        print(f)

        # --- POSITIONS --- #
        if catbool or sim =='eagle':
            pos_list.append(array(f[ptype]['Coordinates'], dtype=float32)) #kpc/h if dalya, mpc if colibre, mpc/h if eagle
        # --- MASSES --- #
        mass_list.append(array(f[ptype]['Masses'], dtype=float32)) #1e10 Msun/h if dalya, Msun if Colibre
        # --- DENSITY --- #
        if sim == 'dalya' or sim =='hello':
            dens_list.append(array(f[ptype]['Density'], dtype=float32))
        elif sim == 'eagle':
            #get density from entropy, need a few params more for eagle to get density from entropy
            aexp_scale_exponent_S  = f[ptype]['Entropy'].attrs['aexp-scale-exponent']
            h_scale_exponent_S     = f[ptype]['Entropy'].attrs['h-scale-exponent']
            CGSConversionFactor_S  = f[ptype]['Entropy'].attrs['CGSConversionFactor']
            aexp_scale_exponent_U  = f[ptype]['InternalEnergy'].attrs['aexp-scale-exponent']
            h_scale_exponent_U     = f[ptype]['InternalEnergy'].attrs['h-scale-exponent']
            CGSConversionFactor_U  = f[ptype]['InternalEnergy'].attrs['CGSConversionFactor']
            aexp_scale_exponent_NH = f[ptype]['ElementAbundance/Hydrogen'].attrs['aexp-scale-exponent']
            h_scale_exponent_NH    = f[ptype]['ElementAbundance/Hydrogen'].attrs['h-scale-exponent']
            CGSConversionFactor_NH = f[ptype]['ElementAbundance/Hydrogen'].attrs['CGSConversionFactor']

            S     = array(f[ptype]['Entropy'], dtype=float32)
            S    *= a**aexp_scale_exponent_S * h**h_scale_exponent_S * CGSConversionFactor_S
            U     = array(f[ptype]['InternalEnergy'], dtype=float32)
            U    *= a**aexp_scale_exponent_U * h**h_scale_exponent_U * CGSConversionFactor_U
            NH    = array(f[ptype]['ElementAbundance/Hydrogen'], dtype=float32)
            NH   *= a**aexp_scale_exponent_NH * h**h_scale_exponent_NH * CGSConversionFactor_NH
            rho_S = ((2./3.) * (U/S))**(3./2.)*NH  #gram / physical cm^-3 
            
            dens_list.append(array(rho_S, dtype=float32))
            del U, S, NH, rho_S

            #dens = concatenate((dens,array(f[ptype]['Density'])))
            #dens = dens*a**(-3.)*h**2.*Msuningram/(kpcincm*1e3)**3./mA
            
            # --- GroupNumber --- #
            if inout and not catbool:
                GN_list.append(array(f[ptype]['GroupNumber']))
        elif sim =='colibre':
            a_exp_D   = float(f[ptype]['Densities'].attrs['a-scale exponent'][0])
            U_L_exp_D = float(f[ptype]['Densities'].attrs['U_L exponent'][0])
            U_M_exp_D = float(f[ptype]['Densities'].attrs['U_M exponent'][0])
            D_factor  = a**a_exp_D * U_L**U_L_exp_D * U_M**U_M_exp_D #a^-3 U_M U_L^-3  [ g pcm^-3 ]

            D     = array(f[ptype]['Densities'], dtype=float32)
            D    *= D_factor

            dens_list.append(D)

            # --- TEMPERATURE --- #
            temp_list.append(array(f[ptype]['Temperatures'], dtype=float32))
            
            if add_props:
                #HI and Z
                Z_list.append(array(f[ptype]['MetalMassFractions'], dtype=float32)/Zsol)
                HMfrac = array(f[ptype]['ElementMassFractions'][:,0], dtype=float32)
                HI_list.append(array(f[ptype]['SpeciesFractions'][:,1])*HMfrac)

        if sim != 'colibre':
            temp_list.append(array(f[ptype]['Temperature'], dtype=float32))

    # --- END OF HDF5 FILES LOOP --- #

    dens = concatenate(dens_list)
    temp = concatenate(temp_list)
    mass = concatenate(mass_list)
    if add_props:
        Z  = concatenate(Z_list)
        HI = concatenate(HI_list)
    if catbool or sim == 'eagle':
        pos = concatenate(pos_list)
    else:
        pos = empty((0,3), dtype=float32)
    if sim == 'eagle' and not catbool and inout:
        GN = concatenate(GN_list)

    # For Eagle, keep only particles in a sub-region 
    if sim == 'colibre':
        pos = pos*1e3*h #kpc/h now
    if sim == 'eagle':
        pos   = pos*1e3 #kpc/h now
        maskE = where((pos[:,0]<=centerH[0]+dx) & (pos[:,0]>=centerH[0]-dx) & (pos[:,1]<=centerH[1]+dx) & (pos[:,1]>=centerH[1]-dx) & (pos[:,2]<=centerH[2]+dx) & (pos[:,2]>=centerH[2]-dx))
        mass  = mass[maskE]
        pos   = pos[maskE]
        dens  = dens[maskE]
        temp  = temp[maskE]
        if inout and not catbool:
            GN    = GN[maskE]

    if verbose:
        np = len(dens)
        print("Number of particles: "+str(np))

    #deal with halo catalog
    if catbool:
        fRv = [0.,1.]
        out     = True
        if not inout:
            exit("--InOut argument mandatory with --Hcat")
        if args.InOut[0] == 'in':
            out = False
            if args.Rvirfrac is not None:
                fRv = [args.Rvirfrac[0],args.Rvirfrac[1]]
            if args.annulus is not None:
                ckpc = array([args.annulus[0],args.annulus[1]])*h #ckpc/h

        if verbose:
            if out:
                print("Plotting for particles outside virial radii")
            else:
                if annulus:
                    print("Plotting for particles within "+str(ckpc[0])+"-"+str(ckpc[1])+" ckpc/h")
                else:
                    print("Plotting for particles in "+str(fRv[0])+"-"+str(fRv[1])+" virial radii")

        if cattype == 'AHF':
            cat = loadtxt(args.Hcat[0],usecols=(5,6,7,11),dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
        elif cattype == 'P2C':
            cat = loadtxt(args.Hcat[0],usecols=(1,2,3,4),delimiter=';',dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
            #convert to ckpc/h
            cat['Hx']   = cat['Hx']*h*1000.
            cat['Hy']   = cat['Hy']*h*1000.
            cat['Hz']   = cat['Hz']*h*1000.
            cat['Rvir'] = cat['Rvir']*h*(1.+z)

        print("Flagging particles...")
    
        Rvir = atleast_1d(cat['Rvir']).astype(float32)
        Hpos = atleast_2d(array([cat['Hx'],cat['Hy'],cat['Hz']], dtype=float32).T)
        if sim == 'eagle':
            maskEH = where((Hpos[:,0]<=centerH[0]+dx) & (Hpos[:,0]>=centerH[0]-dx) & (Hpos[:,1]<=centerH[1]+dx) & (Hpos[:,1]>=centerH[1]-dx) & (Hpos[:,2]<=centerH[2]+dx) & (Hpos[:,2]>=centerH[2]-dx))
            Rvir   = Rvir[maskEH]
            Hpos   = Hpos[maskEH]
        nh   = len(Rvir)
        if verbose:
            print("Kept "+str(nh)+" halos")

        del cat

        tree = cKDTree(Hpos)
        if annulus:
            max_radius = ckpc[1]
        else:
            max_radius = max(fRv[1]*Rvir)

        neighbours = tree.query_ball_point(pos, r=max_radius)
        inorout = zeros(len(pos), dtype=bool) # True if within the requested halo region

        for ip, ih in enumerate(neighbours):
            if len(ih) == 0:
                continue
            ih = array(ih, dtype=int)
            dist = linalg.norm(Hpos[ih]-pos[ip], axis=1)
            if annulus:
                inorout[ip] = any((dist>=ckpc[0]) & (dist<=ckpc[1]))
            else:
                inorout[ip] = any((dist>fRv[0]*Rvir[ih]) & (dist<=fRv[1]*Rvir[ih]))
    
        #restrict to wanted particles
        if out:
            selected = where(inorout==False)
        else:
            selected = where(inorout==True)
        if verbose:
            print("shape(selected) = "+str(shape(selected)))
        dens = dens[selected]  
        temp = temp[selected]
        mass = mass[selected]
        if add_props:
            Z  = Z[selected]
            HI = HI[selected]

    #For eagle, can do in/out Rvir w/o catalog (FRACTION Rvir AND ANNULUS NOT IMPLEMENTED)
    if sim == 'eagle' and not catbool and inout:
        out     = True
        if args.InOut[0] == 'in':
            out = False

        if verbose:
            if out:
                print("Plotting for particles outside virial radii")
            else:
                print("Plotting for particles inside virial radii")

        print("Flagging particles...")

        inorout = where(GN==2**30,0,1) #if GroupNumber = 2^30, particle not in a halo
        
        #restrict to wanted particles
        if out:
            selected = where(inorout==0)
        else:
            selected = where(inorout==1)
        if verbose:
            print("shape(selected) = "+str(shape(selected)))
        dens = dens[selected]
        temp = temp[selected]
        mass = mass[selected]

    #Now let's transform density untis to hydrogen atoms per physical cc
    mH         = 1.6735575e-24 #gram
    mHe        = 6.6464731e-24 #gram
    mA         = mH #0.76*mH+0.24*mHe
    Msuningram = 1.989e43
    kpcincm    = 3.086e21
    if sim == 'dalya' or sim =='hello':
        dens *= a**(-3.)*h**2.*Msuningram/kpcincm**3./mA 
        mass = mass*Msuningram/(h*mA)
    elif sim == 'eagle' or sim == 'colibre':
        dens /= mA
        mass = mass*Msuningram/mA

    #let's get volume as needed for weighted averages
    vol = mass/dens #physical cc
    vol = vol/(kpcincm**3.) #physical kpc to store smaller numbers, and will only be used as weights anyways
        
    #save if wanted
    if savenpy:
        print("Saving density and temperature to "+snappath)
        if inout:
            if out:
                savez(snappath+'/density_temperature_volume_out_Hcat', density=dens, temperature=temp, volume=vol)
            else:
                if catbool:
                    if annulus:
                        if add_props:
                            savez(snappath+'/density_temperature_volume_HI_Z_in_'+str(int(ckpc[0]/h))+'-'+str(int(ckpc[1]/h))+'ckpc_Hcat', density=dens,temperature=temp,volume=vol,HI=HI,metal=Z)
                        else:
                            savez(snappath+'/density_temperature_volume_in_'+str(int(ckpc[0]/h))+'-'+str(int(ckpc[1]/h))+'ckpc_Hcat', density=dens,temperature=temp,volume=vol)

                    else:
                        savez(snappath+'/density_temperature_volume_in_'+str(fRv[0])+'-'+str(fRv[1])+'Rvir_Hcat', density=dens,temperature=temp,volume=vol)
                else:
                    savez(snappath+'/density_temperature_volume_in', density=dens,temperature=temp,volume=vol)
        else:
            savez(snappath+'/density_temperature_volume_all', density=dens,temperature=temp,volume=vol)

# --- End "if --dfile is None"


#If loading precomputed density from file
if args.dfile is not None:
    savenpy = False
    dfile = args.dfile[0]
    data  = load(snappath+'/'+dfile)
    dens  = data['density']
    temp  = data['temperature']
    if add_props:
        HI = data['HI']
        Z  = data['metal']
        HI = where(HI>0.0,HI,1e-6)

    z = 3.5
    print("Redshift = "+str(z))
    if args.annulus is not None:
        ckpc = array([args.annulus[0],args.annulus[1]]) #ckpc in this particular case


#2D HSITOGRAM
if plot and not add_props:
    print("Plotting...")
    binst = linspace(min(log10(temp)),max(log10(temp)),500)
    binsd = linspace(min(log10(dens)),max(log10(dens)),500)
    H, xedges, yedges = histogram2d(log10(dens), log10(temp),bins=[binsd,binst],density=False)
    # some curing of H: remove nan, mask zeros
    H           = array(H)/len(dens)
    #H           = H/amax(H)
    H[isnan(H)] = 0.0
    #H           = H[::-1]
    Hmasked     = ma.masked_where(H==0,H) # Mask pixels with a value of 0

    my_cmap = plt.cm.jet
    plt.pcolormesh(binsd,binst,Hmasked.T,cmap=my_cmap,norm=LogNorm(vmin=1e-6, vmax=1e-3))
    plt.colorbar().set_label(label=r'$N/N_{tot}$',fontsize=14)
    #plt.plot([-1,-1],[2,8],color='grey',ls='--')
    #plt.plot([0,0],[2,8],color='grey',ls='--')
    #plt.plot([1,1],[2,8],color='grey',ls='--')
    plt.ylim(3,5) #plt.ylim(1,8)
    plt.xlim(-3.0,4.0) #-4.5
    fig = plt.gcf()  # get current figure
    fig.set_size_inches(7, 4.5)
    if args.InOut is not None:
        if out:
            plt.text(-4.2,1.5,r'$> R_{\textrm{vir}}$',fontsize=11)
            fout     = snappath+'/TvsD_out.pdf'
        else:
            if args.Rvirfrac is not None:
                fRv = args.Rvirfrac[0]
                plt.text(-4.2,1.5,r'$<$ '+str(fRv)+r'$R_{\textrm{vir}}$',fontsize=11)
            elif annulus and args.dfile is None:
                plt.text(-4.2,1.5,r''+str(round(ckpc[0]/h,0))+r'$\leq r \leq$'+str(round(ckpc[1]/h,0)),fontsize=11)
            elif annulus:
                plt.text(1.5,4.65,r''+str(int(ckpc[0]))+r'$\leq r \leq$'+str(int(ckpc[1])),fontsize=11)
            else:
            	plt.text(-4.2,1.5,r'$\leq R_{vir}$',fontsize=11)
            plt.text(1.5,4.5,r'$12.0 \leq \log M \leq 12.5$',fontsize=11)
            fout     = snappath+'/TvsD_in.pdf'
    else:
        plt.text(-4.2,1.5,r'All particles',fontsize=13)
        fout     = snappath+'/TvsD_all.pdf'
    plt.xlabel(r'$\log(\rho) \;[\textrm{cm}^{-3}]$',fontsize=14, labelpad=12)
    plt.ylabel(r'$\log(T)$ [K]',fontsize=14, labelpad=12)
    #plt.title(r'z = '+str(round(z,2)),fontsize=16)
    plt.text(1.5,4.8,r'z = '+str(round(z,2)),fontsize=11)
    if sim == 'dalya':
        plt.title('DaLya',fontsize=18)
    elif sim == 'hello':
        plt.title('HELLO',fontsize=18)
    elif sim == 'colibre':
        plt.title('COLIBRE m6',fontsize=18)
    plt.plot([-5,5.5],[4,4],'k-',linewidth=0.5)    
    plt.tight_layout()
    plt.savefig(fout, bbox_inches='tight')
    plt.show()

if plot and add_props:
    plt.figure(1)
    # Compute mean Z in each (T, den) bin
    Zmean, binsd, binst, binnumber = binned_statistic_2d(log10(dens), log10(temp), Z,statistic="mean",bins=500)   # "mean", "median", "sum", "std", "count", or custom function
    plt.pcolormesh(binsd, binst, Zmean.T, norm=LogNorm(vmin=1e-3, vmax=1),shading="auto")
    plt.colorbar().set_label(label=r"$Z/Z_\odot$",fontsize=14)
    plt.ylim(3,5)
    plt.xlim(-3.0,4.0)
    fig = plt.gcf()  # get current figure
    fig.set_size_inches(7, 4.5)
    if sim == 'colibre':
        plt.title('COLIBRE m6',fontsize=18)
    plt.xlabel(r'$\log(\rho) \;[\textrm{cm}^{-3}]$',fontsize=14, labelpad=12)
    plt.ylabel(r'$\log(T)$ [K]',fontsize=14, labelpad=12)
    f1 = snappath+'/TvsD_in_Zcolored.pdf'
   # print("Saving figure to "+f1)
   # plt.savefig(f1, bbox_inches='tight')
    
    plt.figure(2)
    # Compute mean HI in each (T, den) bin
    HImean, binsd, binst, binnumber = binned_statistic_2d(log10(dens), log10(temp), HI,statistic="mean",bins=500)   # "mean", "median", "sum", "std", "count", or custom function
    plt.pcolormesh(binsd, binst, HImean.T, norm=LogNorm(vmin=1e-3, vmax=1))#shading="auto")
    plt.colorbar().set_label(label=r"HI mass fraction",fontsize=14)
    plt.ylim(3,5)
    plt.xlim(-3.0,4.0)
    fig = plt.gcf()  # get current figure
    fig.set_size_inches(7, 4.5)
    if sim == 'colibre':
        plt.title('COLIBRE m6',fontsize=18)
    plt.xlabel(r'$\log(\rho) \;[\textrm{cm}^{-3}]$',fontsize=14, labelpad=12)
    plt.ylabel(r'$\log(T)$ [K]',fontsize=14, labelpad=12)
    f2 = snappath+'/TvsD_in_HIcolored.pdf'
    print("Saving figure to "+f2)
    plt.savefig(f2, bbox_inches='tight')

    plt.show()

