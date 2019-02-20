
import nuFATEpy as nuf
import numpy as np
import sys
import time
import math as m

flavor_id = 2
gamma_index = 2.2
h5_filename = "../../resources/NuFATECrossSections.h5"
include_secondaries = False
zenith = 180.* np.pi/180

nufate = nuf.nuFATE(flavor_id, gamma_index, h5_filename, include_secondaries)

print "generate nuFATE"

ntrials = 1000
Na = 6.0221415e23
deltas = np.random.rand(ntrials)
earth_t = nufate.get_earth_column_density(zenith) * Na

print nufate.energy_nodes()
att = nufate.get_relative_attenuation(earth_t)
print att

# do speed test.
start = time.time()
for i in range(ntrials) :
   nufate.set_initial_power_law_flux(2.0 + deltas[i])
   att = nufate.get_relative_attenuation(earth_t)

end = time.time()

print "test1, prosess time was ", end - start

energy_nodes = np.array(nufate.energy_nodes())
unitvec = np.ones(len(energy_nodes))
dummy_initial_flux = np.power(energy_nodes, -2.2*unitvec);

start2 = time.time()
for i in range(ntrials) :
   nufate.set_initial_flux(dummy_initial_flux * np.power(energy_nodes, deltas[i]*unitvec))
   att = nufate.get_relative_attenuation(earth_t)

end2 = time.time();

print "test2, prosess time was ", end2- start2

for i in range(ntrials) :
   nufate.set_initial_power_law_flux(2.0 + deltas[i])
   att = nufate.get_relative_attenuation(earth_t)
   print att

for i in range(ntrials) :
   nufate.set_initial_flux(dummy_initial_flux * np.power(energy_nodes, deltas[i]*unitvec))
   att = nufate.get_relative_attenuation(earth_t)
   print att

