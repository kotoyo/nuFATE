![nuFATE Logo](/resources/nufate.png)

nuFATE is a code that rapidly calculates the attenuated neutrino flux as a function of energy, neutrino flavor and zenith angle, due to the earth's opacity. The software as implemented employs a user-specified power-law isotropic neutrino flux, the STW105 reference earth model, and neutrino-nucleon cross sections computed with the CT10nlo PDF distribution. The attenuation rates can be used to calculate the upgoing nu and nubar fluxes for high-energy neutrino observatories such as IceCube or ANTARES. Full description is available here: https://arxiv.org/abs/1706.09895

Prerequisites
-------------

The following packages are required to use the library, and
can probably be obtained from your favorite package manager:

* numpy: http://www.numpy.org/
* scipy: http://www.scipy.org/
* tables: https://pypi.python.org/pypi/tables

Recommended:
* ipython: http://ipython.org/
* jupyter: http://jupyter.readthedocs.io/
* matplotlib: http://matplotlib.org/

For the C++ version, you will need:

* hdf5 with c bindings: http://www.hdfgroup.org/HDF5/
* gsl (>= 1.15): http://www.gnu.org/software/gsl/
* C++ compiler with C++11 support


Compiling
---------

The Python interface can be installed by simply running:

  python setup.py install

Without write permissions, you can install it using:

  python setup.py install --user

The library can be compiled by running:

	make

An example program demonstrating usage and functionality
can be compiled with the command:

	make examples

The resulting example executables can then be found in the
subdirectory of `examples`

Finally the library can be installed using:

	make install

Compile pybinding (nuFATEpy)
-----------------------------

To compile pybinding, cvmfs service (for icecube) is required.

Also, the main nuFATE project must be installed with 'make install' command in advance.

Load icecube cvmfs environment before typing make command.

        eval `/cvmfs/icecube.opensciencegrid.org/py2-v3/setup.sh`
        cd nuFATE/src/pybinding
        make 

Then, add nuFATE/src/pybinding to your PYTHONPATH to use nuFATEpy.

Example script for using nuFATEpy is nuFATE/src/pybinding/example.py.


How to use nuFATE
--------------------------------------

### 1) Get eigen values.

Following 5 parameters are required to calculate flux after propagation.

        eval : 1D vector, right hand side matrix eigenvalues in unit of cm^2
        evec : 2D vector, right hand side matrix normalized eigenvectors.
        ci   : 1D vector, coordinates of the input spectrum in the eigensystem basis.
        energy_nodes : 1D vector, array containing the energy nodes in GeV.
        phi_0 : 1D vector, input_flux * scaling_flux.

The scaling_flux is by default E^2.0. For most cases the default value works fine. Do not touch it unless you know well about it.

#### For C++

```
#include "nuFATE.h"

....

    // Instanciate object
    int flavor = -2; 
    double gamma = 2.2;
    bool include_secondaries = false;
    std::string file = "../resources/NuFATECrossSections.h5"; 

    //Initialize an instance of the nuFATE class with these parameters.
    nufate::nuFATE object(flavor, gamma, file, include_secondaries);

    //Result is a class that stores the solution to the cascade equation.
    nufate::Result result = object.getEigensystem();
    std::vector<double> eval = result.eval;
    std::shared_ptr<double> evec = result.evec;
    std::vector<double> ci = result.ci;
    std::vector<double> energy_nodes = result.energy_nodes_;
    std::vector<double> phi_0 = result.phi_0_;
```

#### For python

```
    import numpy as np
    import cascade as cas
    import earth

    #Choose the flavor & index you want
    flavor = 2  # 1,2,3 = e, mu, tau; negative sign for antiparticles
    gamma = 2.2  # Power law index of isotropic flux E^-gamma
    ReverseTime = False #You want to go backwards or forward? 
    Efinal = 0.5e9   

    xsecfile = "../../resources/NuFATECrossSections.h5"; 

    #solve the cascade equation once
    eval, evec, ci, energy_nodes, phi_0 = cas.get_eigs(flavor, gamma, xsecfile, ReverseTime, Efinal)
```


### 2) Calculate total number of targets

#### for C++

```
    double Na = 6.0221415e23;
    double zenith = 2.2689280276;
    double t;
    t = object.getEarthColumnDensity(zenith) * Na;
```

#### for python

```
    Na = 6.0221415e23;
    zenith = 2.2689280276;
    t = earth.get_t_earth(zenith) * Na
```

### 3) Calculate flux after propagation.


Basically, all one need to calculate is the following equation.

        phi_out = evec (dot) (ci * exp(eval * t)) 
        
        (dot) : inner product operator
        t : total number of targets


Note that this solution (phi_out) provided by nuFATE is :

        phi_out = scaling_flux * phi. 
        
        where
        phi : flux after propagation
        scaling_flux : E^2 (with default value)

In other words, if you want to get phi, you need to devide the answer (phi_out) with scaling flux.

        phi = phi_out * scaling_flux

Or, if you are interested in ratio of fluxes phi / phi_input, use

        phi_sol = evec (dot) (ci * exp(eval * t)) / phi_0

then the scaling_flux is cancelled out.

In c++ version (and c++ pybinding version) provides a function to get phi_sol : getRelativeAttenuation(t).

For detailed coding, see examples.


### 4) Reverse time propagation.

nuFATE allows to calculate the propagation starting from arrival flux to input flux at Earth's surface.
See example.py in src/python directory for details.


Example
-------

The example script is called example.py. To run it do

python example.py

There is also an iPython notebook (notebook.ipynb) which shows some examples, including plots. This requires the "recommended" packages, above. 

To run the C++ example:

./example.exe


Cautions and some restrictions of NuFATE
-----------------------------------------------

NuFATE provides two modes to calculate attenuated flux. 
The first mode(Mode 1) calculates attenuation of initial neutrino flux which particle type (mentioned as "flavor" in program code) is defined via flavor_id.
The second mode(Mode 2) includes contribution of NuTau's regeneration effect on NuE or NuMu flux, as well as NuTau to NuTau regeneration effect. 
Note that the Mode 2 can be used only when we assume initial flux of NuE:NuEBar:NuMu:NuMuBar:NuTau:NuTauBar as 1:1:1:1:1:1, and if the condition is fulfilled, the authors always recommend using Mode 2 to include NuTau regeneration effect.

The flavor_id must be one of the integer numbers listed below.

* 1 for NuE
* -1 for NuEBar
* 2 for NuMu
* -2 for NuMuBar
* 3 for NuTau
* -3 for NuTauBar

NuEBar, NuTau and NuTauBar may generate neutrinos with other particle types due to Glashow Resonance or Tau regeneration process. Read notes listed below for each case.

### NuEBar

NuEBar may have Glashow Resonance interaction and the outgoing neutrinos from W- decay could be any flavor.

1) NuEBar + e- -> (W-) -> e- + NuEBar
2) NuEBar + e- -> (W-) -> mu- + NuMuBar
3) NuEBar + e- -> (W-) -> tau- + NuTauBar
4) NuEBar + e- -> (W-) -> hadrons

NuFATE takes into account of the first case (NuEBar + e- -> (W-) -> e- + NuEBar) only. In other words, there is no function to get the arrival flux from NuEBar to NuMuBar or NuTauBar. To estimate these contributions, use nuSQuIDS.

Cross sections of Glashow Resonance is hard coded in programs.


### NuTau and NuTauBar

NuTau and NuTauBar may generate tau via CC interaction, and the tau may decay into any flavor:

1) NuTau -> tau -> e + NuEBar + NuTau (18%)
2) NuTau -> tau -> mu + NuMuBar + NuTau (18%)
3) NuTau -> tau -> hadron + NuTau (64%)

With nuFATE Mode 1, the result attenuation ratio includes regeneration from NuTau to NuTau only.

To calculate attenuation ratio of NuE(NuEBar) or NuMu(NuMuBar) with contribution from NuTau's regeneration effect, use cascade_sec.py for python mode or activate "include_secondaries" option in constructor of c++ version. Allowed flavor_ids are -1, 1, -2, and 2.  As mentioned above, with this mode (Mode 2) initial flux is assumed to be 1:1:1:1:1:1 for all particle types.

**Example: Calculate NuEBar's attenuation ratio with NuTau's regeneration effect (Mode 2):**

Set flavor_id as -1, and give initial flux (which is same for all flavors).  
The obtained attenuation ratio has a length of 2N, where N is size of energy_nodes. The first half of the attenuation ratio is for the NuEBar with respect to initial NuEBar flux.  
The second half of attenuation ratio gives NuTauBar's attenuation, and that is same as the attenuation ratio obtained with Mode 1 simulation.


Format of cross sections
------------------------

### Format of NuFATECrossSections.h5

Total cross sections must be charged-current(CC) + neutral-current(NC) cross sections as a function of neutrino energy.
Unit for energy is GeV and for cross sections is cm^2.

Here is the details of each component. () represents shape of matrices, N represents isoscalar particle(0.5(p+n)).

- **total_cross_sections.nuexs(200,)** : NuE-N total cross section (CC + NC)
- **total_cross_sections.nuebarxs(200,)** : NuEBar-N total cross section (CC + NC)
- **total_cross_sections.numuxs(200,)** : NuMu-N total cross section (CC + NC)
- **total_cross_sections.numubarxs(200,)** : NuMuBar-N total cross section (CC + NC)
- **total_cross_sections.nutauxs(200,)** : NuTau-N total cross section (CC + NC)
- **total_cross_sections.nutaubarxs(200,)** : NuTauBar-N total cross section (CC + NC)


- **diffferential_cross_sections.dxsnu (200, 200)** : Nu-N NC differential cross section, the first axis is for primary nu energy and the second axis is for scattered nu energy. Same cross section is used for all flavors.
- **diffferential_cross_sections.dxsnubar (200, 200)** : NuBar-N NC differential cross section, the first axis for primary nubar energy and the second axis for scattered nubar energy. Same cross section is used for all flavors.


- **tau_decay_spectrum.tfull (200,200)** : Differential cross section for NuTau regeneration NuTau -> tau -> NuTau + something. The first axis for primary NuTau energy and the second axis is for regenerated NuTau energy.
- **tau_decay_spectrum.tbarfull (200,200)** : Differential cross section for NuTauBar regeneration NuTauBar -> TauBar -> NuTauBar + something. The first axis is for primary NuTauBar energy and the second axis is for regenerated NuTauBar energy.

** TODO : the following notes may be wrong, confirm it **

- **tau_decay_spectrum_secfull (200,200)** : Differential cross section for NuTau -> tau -> e+NuEBar or mu+NuMuBar. The first axis is for primary NuTau energy and the second axisis for NuEBar or NuMuBar energy.
- **tau_decay_spectrum_secbarfull (200,200)** : Differential cross section for NuTauBar -> taubar -> e+ + NuE or mu+ + NuMu. The first axis is for primary NuTauBar energy and the second axis is for NuE or NuMu energy.

For tau decay component, Appendix A and Table 1 of https://arxiv.org/abs/hep-ph/0005310 is used.


### Format of text-based cross sections
NuFATE accepts text-based cross section files in format of nuSQuIDS cross section.
Currently NuTau regeneration is not supported yet for text-based cross sections. 

Format of total cross section file must have 7 column and N rows where N = number of energy bins.  

- Energy  NuE_Xsec  NuEBar_Xsec  NuMu_Xsec  NuMuBar_Xsec  NuTau_Xsec  NuTauBar_Xsec

Energies must be in GeV and cross sections(Xsec) must be in cm^2.  
Cross section files are separated for CC interaction and NC interaction.


NuFATE uses dsigma/dE differential cross section for NC interaction.  
The differential cross section file must have 8 column and N rows where N = number of energy bins.  

- Energy_in  Energy_out  NuE_Xsec  NuEBar_Xsec  NuMu_Xsec  NuMuBar_Xsec  NuTau_Xsec  NuTauBar_Xsec

Energies must be in GeV and cross sections(Xsec) must be in cm^2.




Citation
--------

If you want cite this work, or want to look at further description
please refer to

High-energy neutrino attenuation in the Earth and its associated uncertainties

Aaron C. Vincent, Carlos A. Arguelles, and A. Kheirandish

arXiv:1706.09895

Contributors
------------

- Aaron C. Vincent
- Carlos A. Arguelles
- Ali Kheirandish
- Ibrahim Safa
- Kotoyo Hoshina

