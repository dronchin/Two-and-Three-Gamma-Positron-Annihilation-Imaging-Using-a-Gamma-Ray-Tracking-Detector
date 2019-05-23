import numpy as np
import matplotlib.pyplot as plt

def cross_section(e1, e2):
    ''' A function that uses the cross section for positronium formation to 
    return the probability that an event with e1 = energy1 and e2 = energy2.
    The energies are scaled between the range [0,1) where 0.5 = 511keV. 
    
    The equation for the cross section comes from Kamińska, et.al. (2016) 
    https://doi.org/10.1140/epjc/s10052-016-4294-3
    '''
    # reject if the energies are 0 because it blows up (rare case in np.random.unifrom())
    if e1==0 or e2==0:
        return 0
    # reject because e1 + e2 > 0.5 has to be true as e3 !> 0.5
    elif e2 < (0.5 - e1) :
        return 0
    # return the cross section of the energies/4 where 4 is the maximum value
    # the cross section outputs with e1,e2 <= 0.5
    else:
        return ((e1 + e2 - 0.5)**2/(e1**2*e2**2))/4

E_all = []
E_onetwo = []
E_last = []
for i in range(1000000): # Note: Does not return the same number of events
    e1, e2 = np.random.uniform(0,0.5,2) # Randomize e1, e2  between [0,0.5)
    prob = cross_section(e1,e2) # Returns the probability based on the cross section
    if np.random.rand() < prob: # Monte Carlo to determine if we accept
        # If accepted, add them to the list of all of the energies
        E_all.append(e1*1022)
        E_all.append(e2*1022)
        E_all.append((1-e1-e2)*1022)
        E_onetwo.append(e1*1022)
        E_onetwo.append(e2*1022)
        E_last.append((1-e1-e2)*1022)
         
# Prove we have the correct distribtuion
plt.hist(E_all,bins=100)
plt.ylabel("counts")
plt.xlabel("Energy (keV)")
plt.title("All energies")
plt.show()

plt.hist(E_onetwo,bins=100, alpha = 0.5)
plt.hist(E_last, bins=100, alpha = 0.5)
plt.ylabel("counts")
plt.xlabel("Energy (keV)")
plt.title("The first two energies seperated from the last")
plt.show()

