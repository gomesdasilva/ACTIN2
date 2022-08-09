.. ACTIN documentation master file, created by
   sphinx-quickstart on Fri Mar 26 10:58:24 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ACTIN: Activity Indices Toolkit
===============================

ACTIN is a python tool that calculates user defined spectroscopic activity indices for different spectrographs.

The method used is based on the original s-index (see `Duncan et al. 1991 <https://ui.adsabs.harvard.edu/abs/1991ApJS...76..383D/abstract>`_), where the flux in the cores of the activity sensitive lines is divided by reference pseudo-continuum regions.
The flux determination, index calculation and errors are decribed in the Appendix A of `Gomes da Silva et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021A%26A...646A..77G/abstract>`_.

Included are some of the most used activity proxies in the optical:

- CaII H&K
- Halpha (using 0.6 and 1.6 angstrom bandpasses)
- NaI D2
- HeI

A description of these indices is provided in `Gomes da Silva et al. (2011) <https://ui.adsabs.harvard.edu/abs/2011A%26A...534A..30G/abstract>`_.
Other indices can be added easily by editing the indices table (see "Adding new indices").

Currently implemented spectrographs are:

- HARPS
- HARPS-N
- ESPRESSO
- SPIRou

New instruments are easy to implement, see "Adding new instruments".

With a simple one line python code, it is possible to run ACTIN in a list of fits files and extract activity time series:

.. code:: python

   from actin import ACTIN
   
   df = ACTIN().run(files, indices=['I_CaII', 'I_NaI'])

where ``files`` is the list of fits file paths to be read, ``ìndices`` is the list of index identification names (as in the indices table), and ``df`` is the pandas dataframe with the output timeseries including the selected headers such as radial-velocity, BJD, and activity data. See the tutorials for more information.

ACTIN is being developed in a `public repository on GitHub`_ so if you have any trouble, open an issue there.

.. _public repository on GitHub: https://github.com/gomesdasilva/ACTIN


Contents
--------

.. toctree::
   :maxdepth: 2

   quickstart
   api
   notebooks/getstarted




License & Attribution
---------------------

Copyright 2018-2021 João Gomes da Silva and contributors.

ACTIN is being developed by João Gomes da Silva in a
`public GitHub repository <https://github.com/gomesdasilva/ACTIN>`_.
The source code is made available under the terms of the MIT license.

If you make use of this code, please cite `Gomes da Silva et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018JOSS....3..667G}>`_:

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
