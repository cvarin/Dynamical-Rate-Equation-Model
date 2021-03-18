"""

"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['ps.useafm'] = True
# plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams.update({'font.size': 14})
import os
import tqdm
import copy

from pyplasma import *
import scipy.constants as c



def fth(material, tau, tolerance=0.01):

    Fmin, Fmax = 0.1, 10 # FIXME: hardcoded limits for minimal and maximal values

    F = 10
    while abs(Fmax-Fmin) > tolerance:

        material_sample = copy.deepcopy(material)

        laser = Laser(wavelength=800*nm, pulse_duration=tau, fluence=F*1e4, t0=0, phase=False)
        dom = Domain()
        dom.add_laser(laser, remove_reflected_part=True)
        dom.add_material(material_sample)
        dom.add_observer(Returner('fi_rate', out_step=1))
        dom.add_observer(Returner('ii_rate', out_step=1))
        dom.add_observer(Returner('rho', out_step=1))
        dom.add_observer(Returner('el_heating_rate', out_step=1))

        def full_ionization(rho):
            return rho > 0.99*material.density
        dom.add_observer(Stopper('rho', full_ionization, verbose=False))

        def is_above_threshold(results):
            absorbed_energy_interband = np.sum(results['fi_rate']*material.bandgap)
            absorbed_energy_intraband = np.sum(results['rho']*results['el_heating_rate']*c.hbar*dom.laser.omega*dom.dt)
            return absorbed_energy_intraband + absorbed_energy_interband > 3e9

        results = dom.run((-2*tau, 2*tau), Nt=500+tau/fs, progress_bar=False)

        if is_above_threshold(results):
            Fmax = F
        else:
            Fmin = F
			
        F = np.exp((np.log(Fmax)+np.log(Fmin))/2)

    return F



if __name__ == '__main__':

	taus = np.logspace(np.log10(20), np.log10(1200), 20)*fs

	materials = {

		'sio2':Material(index=1.45, drude_params={'damping':2e15, 'm_CB':0.75}, 
                        ionization_params={'rate_equation':'dre','bandgap':9*eV,
                                           'density':2.2e28,'cross_section':6.61e-20,
                                           'recombination_rate':1/250e-15}),

		'al2o3':Material(index=1.76, drude_params={'damping':2e15, 'm_CB':0.85}, 
                        ionization_params={'rate_equation':'dre','bandgap':6.5*eV,
                                           'density':2.35e28,'cross_section':1.33e-19}),

		'hfo2':Material(index=2.09, drude_params={'damping':2e15, 'm_CB':0.7}, 
                        ionization_params={'rate_equation':'dre','bandgap':5.1*eV,
                                           'density':2.77e28,'cross_section':1.24e-19}),

		'ta2o5':Material(index=2.1, drude_params={'damping':2e15, 'm_CB':0.85}, 
                        ionization_params={'rate_equation':'dre','bandgap':3.8*eV,
                                           'density':1.12e28,'cross_section':2.50e-19}),
	
		'tio2':Material(index=2.52, drude_params={'damping':2e15, 'm_CB':0.5}, 
                        ionization_params={'rate_equation':'dre','bandgap':3.3*eV,
                                           'density':3.19e28,'cross_section':1.08e-19}),

		# Alternative fit, with recombination
		'tio2_alt':Material(index=2.52, drude_params={'damping':2e15, 'm_CB':0.5}, 
                        ionization_params={'rate_equation':'dre','bandgap':3.3*eV,
                                           'density':3.19e28,'cross_section':1.08e-19,'recombination_rate':1e12}),
	}

	def compute(material,taus):
		results = []
		for tau in tqdm.tqdm(taus, f'Calculating Fth for {material}'):
			results.append(fth(materials[material], tau))
		return results

	# NOTE : Below, we plot the data stored in the fths_... arrays.
    #        If the parameters are changed, these arrays should be removed 
    #        and the computations redone by uncommenting the lines calling
	#        the compute() fonction.

	# fths_sio2 = compute("sio2", taus)
	# print(fths_sio2)
	fths_sio2 = [1.5856061167114663, 1.7583980258575733, 1.9325592266206337, 2.076743273621124, 2.152820461347157, 
	             2.342226586100149, 2.562662455353227, 2.803844384249606, 3.0608344712219995, 3.326386092717056, 
				 3.6149764196179937, 3.902191636323655, 4.193325680003372, 4.468342806945299, 4.70551781378542, 
				 4.916436006163768, 5.159961774877159, 5.458339528628194, 5.85898337414077, 6.3387256192534505]

	# fths_al2o3 = compute("al2o3", taus)
	# print(fths_al2o3)
	fths_al2o3 = [1.1678385645682468, 1.2549684637005685, 1.26062506369433, 1.3607835872433762, 1.4821714294894761, 
				  1.6216642306473654, 1.766323766264441, 1.92388756458243, 2.0861039120183396, 2.2543811464769807, 
				  2.428029205013629, 2.603318734575184, 2.7849937106895744, 2.9659782964256842, 3.1516294950967096, 
				  3.34137932648791, 3.526657476525679, 3.72220922603457, 3.9109760483643834, 4.100086023925162]

	# fths_hfo2 =  compute("hfo2", taus)
	# print(fths_hfo2)
	fths_hfo2 = [0.7136998820779037, 0.7808688975887434, 0.8543594451020262, 0.9263964347542265, 1.0135831333340661, 
				 1.0990453800572955, 1.1890368452040505, 1.289292815223578, 1.3854841522376868, 1.4888521160097425, 
				 1.5927530210435887, 1.696259088072421, 1.8064915626297469, 1.9152548134978824, 2.0214549284871635, 
				 2.1335437974658613, 2.2442654247240537, 2.358080319847817, 2.472102115466655, 2.5800082443957644]

	# fths_ta2o5 = compute("ta2o5", taus)
	# print(fths_ta2o5)
	fths_ta2o5 = [0.4410941012514944, 0.4826071479433916, 0.5280271456481317, 0.5699796823896649, 0.6236226422509303, 
				  0.6762046065333562, 0.726654763727681, 0.7808688975887434, 0.8467093878154393, 0.9017333353397983, 
				  0.969009701214439, 1.0319814123085889, 1.0990453800572955, 1.1625983105011952, 1.2270639228006246, 
				  1.2951041275724522, 1.3669171345425761, 1.4362384940381747, 1.5023039616401745, 1.5714083716097464]

	# fths_tio2 = compute("tio2", taus)
	# print(fths_tio2)
	fths_tio2 = [0.3751619201544639, 0.41046983804365444, 0.4410941012514944, 0.4740031742312117, 0.5093675216780136, 
				 0.5473703262878217, 0.5855690645496172, 0.6236226422509303, 0.6641491558765203, 0.7009759614937345, 
				 0.7465293297167827, 0.7879240830937494, 0.8241677662985539, 0.8698675399968185, 0.9098805364596214, 
				 0.9517340888833216, 0.9955128609158505, 1.0413054106443371, 1.079451477234259, 1.1215139695141063]

	# fths_tio2_alt = compute("tio2_alt", taus)
	# print(fths_tio2_alt)
	fths_tio2_alt = [0.3751619201544639, 0.41046983804365444, 0.4410941012514944, 0.4740031742312117, 0.5093675216780136, 
					 0.5473703262878217, 0.5908597073055873, 0.6292571007878094, 0.6701497732812741, 0.7136998820779037, 
					 0.7600801223642626, 0.8022262617851207, 0.8543594451020262, 0.9017333353397983, 0.9517340888833216, 
					 1.0045073642544629, 1.0602068966819904, 1.1215139695141063, 1.1890368452040505, 1.26062506369433]


	fig = plt.figure(figsize=(6,4))
	colors = plt.cm.terrain(np.linspace(-0.0, 1.0, 19))

	# mero2005 - sio2
	t = [23,30,44,108,150,300,370,635,1100]
	F = [1.7,1.85,2.2,2.8,3.35,4.2,4.45,4.9,6.36]
	plt.scatter(t, F, label=r"$\mathrm{SiO_2}$", marker="o", color=colors[0], s=30)
	# dre
	plt.plot(taus/fs, fths_sio2, color=colors[0])

	# mero2005 - al2o3
	t = [23,30,44,108,150,300,370,635,1100]
	F = [1.32,1.48,1.65,2.08,2.18,2.75,2.86,3.22,4.03]
	plt.scatter(t, F, label=r"$\mathrm{Al_2O_3}$", marker="s", color=colors[2], s=30)
	# dre - al2o3
	plt.plot(taus/fs, fths_al2o3, color=colors[2])

	# #mero2005 - hfo2
	t = [28,32,60,102,150,320,430,740,1200]
	F = [0.85,0.925,1.12,1.31,1.53,1.79,1.9,2.4,2.58]
	plt.scatter(t, F, label=r"$\mathrm{HfO_2}$", marker="D", color=colors[4], s=30)
	# dre - hfo2
	plt.plot(taus/fs, fths_hfo2, color=colors[4])

	#mero2005 - ta2o5
	t = [28,32,60,102,150,320,430,740,1200]
	F = [0.51,0.56,0.64,0.83,0.87,1.06,1.2,1.56,1.81]
	plt.scatter(t, F, label=r"$\mathrm{Ta_2O_5}$", marker="^", color=colors[5], s=30)
	# dre - ta2o5
	plt.plot(taus/fs, fths_ta2o5, color=colors[5])

	#mero2005 - tio2
	t = [28,32,60,102,150,320,430,740,1200]
	F = [0.47,0.49,0.53,0.59,0.63,0.78,0.96,1.05,1.38]
	plt.scatter(t, F, label=r"$\mathrm{TiO_2}$", marker="v", color=colors[7], s=30)
	# dre - tio2
	plt.plot(taus/fs, fths_tio2, color=colors[7])
	plt.plot(taus/fs, fths_tio2_alt,'--', color=colors[7])


	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlim(10, 2000)
	plt.ylim(0.25, 10)
	plt.legend(loc=2, frameon=False, ncol=2, columnspacing=1, handletextpad=0.5)
	plt.xlabel(r"$\tau~\mathrm{[fs]}$")
	plt.ylabel(r"$F_{\mathrm{th}}~\mathrm{[J/cm}^2]$")
	plt.tight_layout()

	plt.savefig("materials_fig.pdf")
	# os.system("pdfcrop materials_fig.pdf materials_fig.pdf > /dev/null")

	# plt.show()