import numpy as np
import matplotlib.pyplot as plt

filename = "out_file.txt"

A = np.loadtxt(filename)

#want to pick out every 10th index

t = A[:,0]
E_perspin_perspin = A[:,1]
E_var_pertemp_pertemp = A[:,2]
M_var_pertemp = A[:,3]
M_abs_perspin_perspin = A[:,4]
accept_config = A[:,5]
time = np.linspace(0, len(E_perspin_perspin)*100, len(E_perspin_perspin), endpoint = True)

#program condition to separate semilogxs based on value of temperature, yo 

E_list = []

for E_val in E_perspin_perspin[1e4:]:
	if not (E_val in E_list):
		E_list.append(E_val)


"""
f, axarr = plt.subplots(4, sharex = True)
axarr[0].semilogx(time, E_perspin_perspin, label = r"$\frac{\langle E \rangle}{N^2}$")
axarr[1].semilogx(time, E_var_pertemp_pertemp, label = r"$\frac{C_v}{T^2}$")
axarr[2].semilogx(time, M_var_pertemp, label = r"$\frac{\chi}{T}$")
axarr[3].semilogx(time, M_abs_perspin_perspin, label = r"$\frac{\langle |M| \rangle }{N^2}$")

f.subplots_adjust(hspace=0.1)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

axarr[0].set_title("semilogx of various expectation-values and means")

[axarr[x].legend(loc = "best", prop={'size':24}) for x in range(4)]
plt.show()


plt.semilogx(time, accept_config)
plt.xlabel("Time- monte carlo cylces")
plt.title("Number of accepted configurations")
plt.show()

"""