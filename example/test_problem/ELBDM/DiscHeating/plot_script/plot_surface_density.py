import argparse
import sys
import yt
import matplotlib.pyplot as plt
import math
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the surface density of the particles' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

center_mode = 'max'
dpi         = 150
field       = 'ParDens'
width = 1.5
disc_thickness = 0.1
Disc_Decay_R = 0.5
TotM = 6.594E-9 * (2E6)

r_max = 1.2
Const = TotM/(Disc_Decay_R-(Disc_Decay_R+r_max)*math.exp(-r_max/Disc_Decay_R))
Const /= (2*3.141*Disc_Decay_R)
yt.enable_parallelism()
ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   my_disk = ds.disk([1.5, 1.5, 1.5], [0.,0.,1.],(1.3, 'code_length'),(disc_thickness, 'code_length') )
   prof = yt.ProfilePlot( my_disk, ("index","cylindrical_r"), 'ParDens', weight_field='cell_volume', n_bins=64, x_log = False )
   surface_dens = prof.profiles[0]["ParDens"].in_units("code_mass/(code_length*code_length*code_length)").d
   surface_dens *= disc_thickness
   logSurDens = np.log10(surface_dens)
   radius = prof.profiles[0].x.in_units("code_length").d
   sumM = 2*3.1416*radius*1.3/64*surface_dens
   print('TotM=',np.sum(sumM))
#   print('Sur=',surface_dens)
#   print('log=',logSurDens)
#   coefficient = np.polyfit(radius[0:-4],logSurDens[0:-4],1)
#   print(coefficient)
   
#   print(my_disk.quantities.total_mass())
   theoretical_sd = np.zeros(len(surface_dens))
   for i in range(len(surface_dens)):
      theoretical_sd[i] = Const*math.exp(-radius[i]/Disc_Decay_R)
   ## plot
   plt.semilogy(radius,surface_dens,'ro')
   plt.plot(radius, theoretical_sd, 'b--')
   plt.ylabel('surface density(code_mass/code_length^2)')
   plt.xlabel('radius(code_length)')
   plt.ylim((10E-6,10E-2))
   plt.title('Surface Density Plot'+str(ds))
   dataname = 'Surface_Density_'+str(ds)+'.png'
   plt.show()
#   plt.savefig(str(dataname))


#  f, ax = plt.subplots(2,1)
#  f.subplots_adjust( wspace = 0.4)
#  [f.axes[t].set_xlim(0.0, 14.0) for t in range(0,4,1)]
#  ax[1][0].set_xlabel( "$\mathrm{Cylindrical\ radius\ [code_lenght]}$", fontsize='large' )
#  ax[0][0].plot( radius, surface_dens, 'r-o', lw=2, mec='none', ms=markersize )
#  ax[0][0].set_yscale( 'log', nonposy='clip' )
#   ax[0][0].set_ylim( 1.0e0, 2.0e3 )
#  ax[0][0].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )
#  ax[0][0].set_ylabel( "$\mathrm{\Sigma_{gas}\ [M_{\odot}/pc^2]}$", fontsize='large' )

#  prof.set_unit( field, 'code_mass/(code_length*code_length)' )
#  prof.set_ylim( field, 1.0e-6, 1.0e0 )
#  prof.save( mpl_kwargs={"dpi":dpi} )

#  add title
#  time = ds.current_time.in_units('Myr')
#  plt.suptitle( "t = %6.2f %s"%(time.d, time.units), fontsize='large' )

#  show/save figure
#  plt.savefig( fileout+'_'+ds.basename+".png", bbox_inches='tight', pad_inches=0.05, dpi=dpi )
#  plt.show()

