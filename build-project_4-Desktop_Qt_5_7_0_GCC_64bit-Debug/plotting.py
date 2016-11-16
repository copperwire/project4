import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

filename = "out_1e6_T_24_order.txt"

A = np.loadtxt(filename)

#want to pick out every 10th index

E_perspin_perspin = A[:,1]
E_var_pertemp_pertemp = A[:,2]
M_var_pertemp = A[:,3]
M_abs_perspin_perspin = A[:,4]
accept_config = A[:,5]
temp = A[:,0]
time = np.linspace(0, len(E_perspin_perspin)*int(1e6/float(len(temp))), len(E_perspin_perspin), endpoint = True)

#program condition to separate semilogxs based on value of temperature, yo 

def equal(a, b, eps = 1e-6):
	if  abs(a - b) < eps:
		return True
	else:
		return False

E_list = [E_perspin_perspin[int(1e6/float(len(temp)))]]
list_vals = [1]


for E_val in E_perspin_perspin[int(1e6/float(len(temp))) + 1:]:
	test_list = [equal(E, E_val) for E in E_list]

	if not (True in test_list):
		E_list.append(E_val)
		list_vals.append(1)
	else:
		list_vals = [iter_val + bol_val * 1 for bol_val, iter_val in zip(test_list, list_vals)] 

list_vals = np.array([x for y, x in sorted(zip(E_list, list_vals))])

s = sum(list_vals)

norm = [float(i)/s for i in list_vals]

m = max(list_vals)

#norm = [float(i)/m for i in list_vals]

E_list = np.array(sorted(E_list))

print(sum(norm))
print(np.average(E_list))
print(np.var(E_list))

plt.plot(E_list, norm)
plt.title("Energy distribution for a temperature of T = 2.4")
plt.xlabel(r"$E$")
plt.ylabel("P(E)")
plt.show()


"""
plt.semilogy(temp, accept_config)
plt.xlabel("T")
plt.title(" Total accepted configurations as a function of temperature")
plt.show()
"""

"""
f, axarr = plt.subplots(2, 2)

i = 0
j = 0
for word in ["random", "order"]:
	for temp_str in ["1", "24"]:
		filename = "out_1e6_T_"+temp_str+"_"+word+".txt"
		A = np.loadtxt(filename)

		#want to pick out every 10th index

		E_perspin_perspin = A[:,1]
		E_var_pertemp_pertemp = A[:,2]
		M_var_pertemp = A[:,3]
		M_abs_perspin_perspin = A[:,4]
		accept_config = A[:,5]
		temp = A[:,0]
		time = np.linspace(0, len(E_perspin_perspin)*int(1e6/float(len(temp))), len(E_perspin_perspin), endpoint = True)

		axarr[i, j].plot(time, accept_config)
		if temp_str == "24":
			temp_str = "2.4"
		axarr[i, j].set_title("Temperature = {}, Initial configuration = {}".format(temp_str, word) )
		axarr[i, j].set_xlabel("Time")
		i += 1
	j+=1
	i = 0

plt.show()
"""
#print(sum(E_var_pertemp_pertemp)/float(len(E_var_pertemp_pertemp	)))


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


"""



"""
plt.semilogx(time, accept_config)
plt.xlabel("Time- monte carlo cylces")
plt.title("Number of accepted configurations")
plt.show()
"""

