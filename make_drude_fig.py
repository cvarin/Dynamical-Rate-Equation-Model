
import matplotlib.pyplot as plt
plt.rcParams['ps.useafm'] = True
# plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams.update({'font.size': 14})
import os

from pyplasma import *


if __name__ == '__main__':

    f = plt.figure(figsize=(6,3))
    colors = plt.cm.terrain(np.linspace(0.0, 1.0, 24))
    ax = f.add_subplot(1,1,1)
   
    gammas = np.logspace(-2,2,1000)
    omega = 1

    # Scaling of E_p vs gamma
    def f2(gamma):
        return 1.0/(gamma**2.0+omega**2.0)
    # Scaling of gamma_lh vs gamma
    def f1(gamma):
        return 2.0*gamma*f2(gamma)

    ax.semilogx(gammas,f2(gammas),color="darkred",lw=2,label=r"$\mathcal{E}_p$")
    ax.semilogx(gammas,f1(gammas),color=colors[2],lw=2,label=r"$\gamma_\mathrm{lh}$")
    
    ax.text(2,0.25,r"$\mathcal{E}_p$",size=16)
    ax.text(2,0.85,r"$\gamma_\mathrm{lh}$",size=18)
    ax.set_xlabel(r"$\gamma/\omega$")
    ax.set_ylabel(r"Normalized amplitude")
    
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels([r"$0$",r"$0.5$",r"$1$"])
    
    plt.tight_layout()
    plt.savefig("drude_fig.pdf")
    # os.system("pdfcrop drude_fig.pdf drude_fig.pdf > /dev/null")

    # plt.show()
