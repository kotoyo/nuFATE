#===================================================
# test 1 : load test
# Is nuFATE object loading cross sections correctly?
#===================================================

import nuFATEpy as nuf
import numpy as np
import sys
import tables
import time

flavor_id = 2 # NuMu
gamma_index = 2.2
include_secondaries = False
zenith = 180.* np.pi/180
Na = 6.0221415e23

h5_filename = "../NuFATECrossSections.h5"

nufate = nuf.nuFATE(flavor_id, gamma_index, h5_filename, include_secondaries)
earth_t = nufate.get_earth_column_density(zenith)

ntrials = 100
np.random.seed(36)
delta_gammas = np.random.random(ntrials)


#---------------------------
# trial 1 
# call constructor everytime
#---------------------------

start_time1 = time.time()

for i in range(ntrials) :
    nufate = nuf.nuFATE(flavor_id, gamma_index+delta_gammas[i], h5_filename, include_secondaries)
    nufate.get_relative_attenuation(earth_t * Na)

end_time1 = time.time()

print "bruteforce repeats : computing time %f sec" % (end_time1 - start_time1) 


#---------------------------
# trial 2 
# set_initial_power_low_flux
#---------------------------

nufate = nuf.nuFATE(flavor_id, gamma_index+delta_gammas[i], h5_filename, include_secondaries)

start_time2 = time.time()

for i in range(ntrials) :
    nufate.set_initial_power_law_flux(gamma_index + delta_gammas[i])
    nufate.get_relative_attenuation(earth_t *Na)

end_time2 = time.time()

print "set_initial_power_law_flux repeats : computing time %f sec" % (end_time2 - start_time2) 

#---------------------------
# trial 3 
# set_initial_flux
#---------------------------

energy_nodes = nufate.energy_nodes()

initial_flux = []
for i in range(ntrials):
   initial_f = [e**-(gamma_index + delta_gammas[i]) for e in energy_nodes]
   initial_flux.append(initial_f)

start_time3 = time.time()

for i in range(ntrials) :
    initial_f = [e**-(gamma_index + delta_gammas[i]) for e in energy_nodes]
    nufate.set_initial_flux(initial_f)
    nufate.get_relative_attenuation(earth_t *Na)

end_time3 = time.time()

print "set_initial_flux repeats : computing time %f sec" % (end_time3 - start_time3) 

#----------------------------
# tests. 
# all three methods should give
# the same answer.
#----------------------------

atten_1 = []
for i in range(ntrials) :
    nufate = nuf.nuFATE(flavor_id, gamma_index+delta_gammas[i], h5_filename, include_secondaries)
    atten_1.append(nufate.get_relative_attenuation(earth_t * Na))

nufate = nuf.nuFATE(flavor_id, gamma_index+delta_gammas[i], h5_filename, include_secondaries)

atten_2 = []
for i in range(ntrials) :
    nufate.set_initial_power_law_flux(gamma_index + delta_gammas[i])
    atten_2.append(nufate.get_relative_attenuation(earth_t *Na))

atten_3 = []
for i in range(ntrials) :
    initial_f = [e**-(gamma_index + delta_gammas[i]) for e in energy_nodes]
    nufate.set_initial_flux(initial_f)
    atten_3.append(nufate.get_relative_attenuation(earth_t *Na))

atten_1 = np.array(atten_1)
atten_2 = np.array(atten_2)
np.testing.assert_allclose(atten_2, atten_1, rtol=1e-5, atol=0, err_msg="test1, bruteforce and set_initial_power_low_flux didn't match")
atten_3 = np.array(atten_3)
np.testing.assert_allclose(atten_3, atten_1, rtol=1e-5, atol=0, err_msg="test2, bruteforce and set_initial_flux didn't match")

print "all tests passed!"
