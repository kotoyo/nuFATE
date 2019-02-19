
import nuFATEpy as nuf
import numpy as np
import sys

flavor_id = 2
gamma_index = 2.2
h5_filename = "../../resources/NuFATECrossSections.h5"
include_secondaries = False
zenith = 180.* np.pi/180

nufate = nuf.nuFATE(flavor_id, gamma_index, h5_filename, include_secondaries)

print "generate nuFATE"

# do speed test.

ntrials = 10000
deltas = np.random.rand(ntrials)
earth_t = nufate.get_earth_column_density(zenith)

for i in range(ntrials) :
   nufate.set_initial_power_law_flux(2.0 + deltas[i]);
   att = nufate.get_relative_attenuation(earth_t)
   

print "done!"



