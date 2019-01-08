![nuFATE Logo](/resources/nufate.png)

nuFATE is a code that rapidly calculates the attenuated neutrino flux as a function of energy, neutrino flavor and zenith angle, due to the earth's opacity, for neutrino energies above 1 TeV. The software as implemented employs a user-specified power-law isotropic neutrino flux, the STW105 reference earth model, and neutrino-nucleon cross sections computed with the CT10nlo PDF distribution. The attenuation rates can be used to calculate the upgoing nu and nubar fluxes for high-energy neutrino observatories such as IceCube or ANTARES. Full description is available here: https://arxiv.org/abs/1706.09895

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


How to use nuFATE
--------------------------------------

## 1) Get eigen values.

Following 5 parameters are required to calculate flux after propagation.

        eval : 1D vector, right hand side matrix eigenvalues in unit of cm^2
        evec : 2D vector, right hand side matrix normalized eigenvectors.
        ci   : 1D vector, coordinates of the input spectrum in the eigensystem basis.
        energy_nodes : 1D vector, array containing the energy nodes in GeV.
        phi_0 : 1D vector, input_spectrum * pedestal_spectrum.

The pedestal_spectrum must be chosen to make phi_0 to be flat. 

For example, if input_spectrum is E^-2.2 (gamma is 2.2), pedestal_spectrum may be E^2.0 (pedestal_index = 2.0).

### For C++

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

### For python

```
    import numpy as np
    import cascade as cas
    import earth

    #Choose the flavor & index you want
    flavor = 2  # 1,2,3 = e, mu, tau; negative sign for antiparticles
    gamma = 1.2  # Power law index of isotropic flux E^-gamma
    ReverseTime = False #You want to go backwards or forward? 
    Efinal = 0.5e9   

    xsecfile = "../../resources/NuFATECrossSections.h5"; 

    #solve the cascade equation once
    eval, evec, ci, energy_nodes, phi_0 = cas.get_eigs(flavor, gamma, xsecfile, ReverseTime, Efinal)
```


## 2) Calculate total number of targets

### for C++

```
    double Na = 6.0221415e23;
    double zenith = 2.2689280276;
    double t;
    t = object.getEarthColumnDensity(zenith) * Na;
```

### for python

```
    Na = 6.0221415e23;
    zenith = 2.2689280276;
    t = earth.get_t_earth(zenith) * Na
```

## 3) Calculate flux after propagation.


Basically, all one need to calculate is the following equation.

        phi_out = evec (dot) (ci * exp(eval * t)) 
        
        (dot) : inner product operator
        t : total number of targets


Note that this solution (phi_out) provided by nuFATE is :

        phi_out = E^(pedestal_index) * phi. 
        
        where
        phi : flux after propagation

In other words, if you want to get phi, you need to devide the answer (phi_out) with E^(pedestal_index).

        phi = phi_out * E^(-pedestal_index)

Or, if you are interested in ratio of fluxes phi / phi_input, use

        phi_sol = evec (dot) (ci * exp(eval * t)) / phi_0

then the pedestal_flux is cancelled out.

For detailed coding, see examples.


## 4) Reverse time propagation.

nuFATE allows to calculate the propagation starting from arrival flux to input flux at Earth's surface.
See example.py in src/python directory for details.


Example
-------

The example script is called example.py. To run it do

python example.py

There is also an iPython notebook (notebook.ipynb) which shows some examples, including plots. This requires the "recommended" packages, above. 

To run the C++ example:

./example.exe

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

