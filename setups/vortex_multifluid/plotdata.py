from pylab import *


### Data
path    = "../../outputs/vortex_multifluid/"
nx      = 128
ny      = 128
n       = 400
rhod1   = fromfile( path+"dust1dens{:}.dat".format(n)).reshape(ny,nx)
rhod2   = fromfile( path+"dust2dens{:}.dat".format(n)).reshape(ny,nx)
rhog    = fromfile( path+"gasdens{:d}.dat".format(n)).reshape(ny,nx)
rhog0   = fromfile( path+"gasdens0.dat").reshape(ny,nx)

print(rhod1.max(), rhod1.min())
print(rhod2.max(), rhod2.min())
print(rhog.max(), rhog.min())
## Plot     
fig = figure(figsize=(14,5.5))

ax1 = fig.add_axes([0.05,0.08,0.275,0.7])
ax2 = fig.add_axes([0.05,0.9 ,0.275,0.02])

ax3 = fig.add_axes([0.4 ,0.08, 0.275,0.7])
ax4 = fig.add_axes([0.4,0.9,0.275,0.02])

ax5 = fig.add_axes([0.75,0.08,0.275,0.7])
ax6 = fig.add_axes([0.75,0.9,0.275,0.02])

vmax    = 1
vmin    = -3
cmap    = 'Spectral'

im1 = ax1.imshow(log10(rhod1), origin='lower', aspect='auto', cmap=cmap,vmax=log10(rhod1).max(),vmin=log10(rhod1).min())
cb1  = colorbar(im1, cax=ax2,orientation="horizontal")
cb1.ax.set_xlabel(r"$\log_{10}(\rho_{\rm d1}/\rho_{\rm g})$", fontsize=14)

im2 = ax3.imshow(log10(rhod2), origin='lower', aspect='auto', cmap=cmap, vmax=log10(rhod2).max(),vmin=log10(rhod2).min())
cb2  = colorbar(im2, cax=ax4,orientation="horizontal")
cb2.ax.set_xlabel(r"$\log_{10}(\rho_{\rm d2}/\rho_{\rm g})$", fontsize=14)

im3 = ax5.imshow(log10(rhog/rhog0), origin='lower', aspect='auto', cmap=cmap, vmax=log10(rhog/rhog0).max(), vmin=log10(rhog/rhog0).min())
cb3  = colorbar(im3, cax=ax6,orientation="horizontal")
cb3.ax.set_xlabel(r"$\log_{10}(\rho_{\rm g}/\rho_{\rm g0})$", fontsize=14)

#tight_layout()
savefig("plotdata.png",  bbox_inches='tight', dpi=250)
