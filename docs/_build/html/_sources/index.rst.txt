ACTIN: Activity Indices Toolkit
===============================

ACTIN is a python tool that calculates user defined spectroscopic activity indices for different spectrographs.

The method used is based on the original s-index (see `Duncan et al. 1991 <https://ui.adsabs.harvard.edu/abs/1991ApJS...76..383D/abstract>`_), where the flux in the cores of the activity sensitive lines is divided by reference pseudo-continuum regions.
The flux determination, index calculation and errors are described in the Appendix A of `Gomes da Silva et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021A%26A...646A..77G/abstract>`_.
The errors are the photon noise, however there is an option to include extra noise, which is added in quadrature to the overall errors.

Included are some of the most used activity proxies in the optical:

- CaII H&K (S-index)
- H$\alpha$ (using 0.6 and 1.6 angstrom bandpasses)
- NaI D2
- HeI

A description of these indices is provided in `Gomes da Silva et al. (2011) <https://ui.adsabs.harvard.edu/abs/2011A%26A...534A..30G/abstract>`_.
Other indices can be added easily by editing the indices table (see `The indices table <file:///Users/jgsilva/Astrophysics/Packages/ACTIN2/docs/_build/html/add_index.html>`_).

Implemented spectrographs:

- HARPS
- HARPS-N
- ESPRESSO

New instruments are easy to implement (see `Adding new spectrographs <file:///Users/jgsilva/Astrophysics/Packages/ACTIN2/docs/_build/html/add_spec.html>`_).
Activity indices can also be calculated for any spectrum if the wavelength and flux are available (see `Using ACTIN with any spectra <file:///Users/jgsilva/Astrophysics/Packages/ACTIN2/docs/_build/html/calc_act_general.html>`_).

With a simple one line python code, it is possible to run ACTIN on a list of fits files and extract a list of activity indices:

.. code:: python

   from actin2 import ACTIN
   actin = ACTIN()
   
   df = actin.run(files, indices=['I_CaII', 'I_NaI'])

where ``files`` is the list of fits file paths to be read, ``ìndices`` is the list of index identification names (as in the indices table), and ``df`` is the pandas DataFrame with the output time series including the selected headers such as radial-velocity, Julian date, and activity data. See the tutorials for more information.

ACTIN is being developed in a `public repository on GitHub`_, if you have any trouble, please open an issue there.

.. _public repository on GitHub: https://github.com/gomesdasilva/ACTIN


.. _installation:

Installation
------------

Clone the `github repository <https://github.com/gomesdasilva/ACTIN>`_ to a directory of your choice and install via:

.. code:: bash

   cd path/to/directory
   python setup.py install


Contents
--------

.. toctree::
   :maxdepth: 2

   getstarted
   read_plot_spectrum
   calc_act_general
   CaII_to_rhk
   add_index
   mod_indices_table
   add_spec
   api



License & Attribution
---------------------

ACTIN is being developed in a
`public GitHub repository <https://github.com/gomesdasilva/ACTIN2>`_.
The source code is made available under the terms of the MIT license.

If you make use of this code, please cite `Gomes da Silva et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018JOSS....3..667G}>`_ and `Gomes da Silva et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021A%26A...646A..77G/abstract>`_:

.. code-block:: tex

    @ARTICLE{2018JOSS....3..667G,
       author = {{Gomes da Silva}, Jo{\~a}o and {Figueira}, Pedro and {Santos}, Nuno and {Faria}, Jo{\~a}o},
        title = "{ACTIN: A tool to calculate stellar activity indices}",
      journal = {The Journal of Open Source Software},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2018,
        month = nov,
       volume = {3},
       number = {31},
        pages = {667},
          doi = {10.21105/joss.00667},
   archivePrefix = {arXiv},
       eprint = {1811.11172},
   primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2018JOSS....3..667G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }



   @ARTICLE{2021A&A...646A..77G,
       author = {{Gomes da Silva}, J. and {Santos}, N.~C. and {Adibekyan}, V. and {Sousa}, S.~G. and {Campante}, T.~L. and {Figueira}, P. and {Bossini}, D. and {Delgado-Mena}, E. and {Monteiro}, M.~J.~P.~F.~G. and {de Laverny}, P. and {Recio-Blanco}, A. and {Lovis}, C.},
        title = "{Stellar chromospheric activity of 1674 FGK stars from the AMBRE-HARPS sample. I. A catalogue of homogeneous chromospheric activity}",
      journal = {\aap},
     keywords = {catalogs, stars: activity, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2021,
        month = feb,
       volume = {646},
          eid = {A77},
        pages = {A77},
          doi = {10.1051/0004-6361/202039765},
   archivePrefix = {arXiv},
       eprint = {2012.10199},
   primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021A&A...646A..77G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }