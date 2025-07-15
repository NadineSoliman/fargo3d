from pylab import *
from reader import Fields

path = "../../outputs/soundwave/"

def eigen(fr, fi,  kmode,  x, lambda_n, time):
    li = lambda_n.imag
    lr = lambda_n.real

    sol = exp(-lr*time) * ( fr*cos(kmode*x - li*time) - fi*sin(kmode*x - li*time) )
    return sol

dt        = 0.1
ntot      = 100
amplitude = 1e-4
NFLUIDS   = 5

Ts        = logspace(-1,0,4)
print("Ts=", Ts)
eps       = linspace(0.1,0.5,4)
L         = 1.0
kmode     = 2.*pi/L
rhog      = 1.0
delta_rho = 1.0
lambda_n  = 0.9124135035415595-5.49379966793634j

vgas      = -1j*lambda_n/kmode *delta_rho/rhog
vdust     = -1j*lambda_n/( kmode*(1-lambda_n*Ts) ) *delta_rho/rhog
rdust     = eps*delta_rho/(1.-lambda_n*Ts)
rgas      = 1.0

time    = []
F3Dvgas = []
F3Drgas = []


x       = loadtxt(path+"domain_z.dat")[3:-3]
nx      = 0


amp  = 1e-4

time     = [] 
for i in range(ntot):
    F3Drgas.append(fromfile(path+"gasdens{0:d}.dat".format(i))[nx]-1.0)
    F3Dvgas.append(fromfile(path+"gasvz{0:d}.dat".format(i))[nx])
    time.append(i*dt)

time  = array(time)
time2 = linspace(time.min(),time.max(),1000)


plot(time2, eigen( vgas.real , vgas.imag , kmode, x[nx] , lambda_n, time2), 'k-')
plot(time, array(F3Dvgas)/amp,'o')

for m in range(NFLUIDS-1):
    F3Dvdust = []
    F3Drdust = []
    for i in range(ntot):
        print(path+"dust{0:d}vz{1:d}.dat".format(m+1,i))
        F3Dvdust.append(fromfile(path+"dust{0:d}vz{1:d}.dat".format(m+1,i))[nx])
        F3Drdust.append(fromfile(path+"dust{0:d}dens{1:d}.dat".format(m+1,i))[nx]-eps[m])
    plot(time2, 5.*eigen( vdust[m].real , vdust[m].imag , kmode, x[nx] , lambda_n, time2), 'k-')
    plot(time, array(F3Dvdust)/amp*5,'o', label=r"$ \epsilon = {:2.2f} - T_s = {:2.2f}$".format(eps[m],Ts[m]))
legend()

xlim(0,2.5)
show()
