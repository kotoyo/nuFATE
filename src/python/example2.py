""" Brief example of how to use nuFATE
"""

import numpy as np
import cascade as cas
import earth
import sys

#Choose the flavor & index you want
flavor = 2  # 1,2,3 = e, mu, tau; negative sign for antiparticles
#gamma = 2.0  # Power law index of isotropic flux E^-gamma
gamma = "../../resources/1PeV_flux.txt"  # Power law index of isotropic flux E^-gamma
ReverseTime = False #You want to go backwards or forward? True for backwards, false for forward in time
Efinal = 0.5e9 #If you're going backwards in time, set the final neutrino energy. The solution in this case returns a pdf
               # of neutrino energies that would give you Efinal, if propagated forwards. 

#gamma = 'data/phiHGextrap.dat' #This is an example Honda Gaisser atmospheric flux. You can use this or add your own file, being careful to follow the energy spacing

#solve the cascade equation once
w, v, ci, energy_nodes, phi_0 = cas.get_eigs(flavor, gamma, "../../resources/NuFATECrossSections.h5", ReverseTime, Efinal)


#this function just interpolates the solution
def get_att_value(w, v, ci, energy_nodes, zenith, E):
    Na = 6.0221415e23
    logE = np.log10(E)
    t = earth.get_t_earth(zenith) * Na
    print "t is ", t
    #t = 6.59272246299e+33 * Na

    # g/ cm^2
    #    phi = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(-2) #this is the attenuated flux
    if(ReverseTime):
        t = -1.*t
        phiend = np.dot(v, (ci * np.exp(w * t))) * energy_nodes**(-2)
        #print phisol
    else:
        phiend = np.dot(v, (ci * np.exp(w * t))) * energy_nodes**(-2) # this is arrival energies
    return np.interp(logE, np.log10(energy_nodes), phiend), phiend


#specify a zenith angle and energy you are interested in
#zenith = np.radians(120.) # zenith angle in radians
zenith = 3.1415# zenith angle in radians
E = 100e3  #GeV
att,phiend = get_att_value(w, v, ci, energy_nodes, zenith, E)

print "Flux at E  =", E, " GeV , zenith = ", np.degrees(zenith), " degrees will be attenuated by a factor of ", att

data = ""
for i, e in enumerate(energy_nodes) :
    data += "%f, %e\n" % (np.log10(e), phiend[i])

f = open("example_out2_python.txt","w")
f.write(data)
f.close()

sys.exit(0)

#done

# The next section shows how to include secondary electron, muon neutrinos
import cascade_secs as csx
w, v, ci, energy_nodes, phi_0 = csx.get_eigs(
    flavor, gamma, "../../resources/NuFATECrossSections.h5")


def get_att_value_secs(w, v, ci, energy_nodes, zenith, E):
    Na = 6.0221415e23
    logE = np.log10(E)
    t = earth.get_t_earth(zenith) * Na
    # g/ cm^2
    #    phi = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(-2) #this is the attenuated flux
    phisol = np.dot(v, (ci * np.exp(w * t))
                   ) / phi_0  #this is phi/phi_0, i.e. the relative attenuation
    phisol = phisol[0:200]  #the non-tau bit.
    return np.interp(logE, np.log10(energy_nodes), phisol)


att = get_att_value_secs(w, v, ci, energy_nodes, zenith, E)

print "If I include secondaries, that factor becomes: ", att

#done done
