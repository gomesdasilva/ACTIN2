.. _quickstart:


Getting started
---------------

.. _installation:

Installation
++++++++++++

Clone the `github repository <https://github.com/gomesdasilva/ACTIN>`_ to a directory of your choice and install via:

.. code:: bash

   cd path/to/directory
   python setup.py install


Example: Extracting a spectrum
++++++++++++++++++++++++++++++


In the following code we will extract the spectrum of the fits ``file`` and plot it:

.. code:: python

    import os
    import matplotlib.pylab as plt
    from actin import ACTIN

    file = os.path.join(actin.dir, "test", "HARPS.2003-02-18T08:28:32.570_s1d_A.fits"))

    actin = ACTIN()

    spec = actin.ReadSpec(file)

    wave = spec.spectrum['wave']
    flux = spec.spectrum['flux']

    plt.plot(wave, flux)
    plt.show()

However, ``ACTIN`` already comes with a handy built-in plotting tool for spectra, so there is no need to extract the wavelength and flux values explicitly. We could just have written:

.. code:: python

    spec.plot(show=True)
    


.. module:: actin2



