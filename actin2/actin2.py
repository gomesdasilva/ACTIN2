import sys, os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import tqdm
import time

from .spectrographs._spec_tools import printif
from .calcindices import CalcIndices
from .readspec import ReadSpec
from .indtable import IndTable


# TODO: add option/class/function to plot lines from IntTable (with the bandpass)

# TODO POST-PROCESSING: option to select HARPS-pre/HARPS-pos (option)
# TODO POST-PROCESSING: add correction for CCF FWHM and contrast drift (option)?

class ACTIN:
    """The ACTIN class. Reads fits files and calculates activity indices.

    Attributes:
        ReadSpec (actin2.ReadSpec) : Object that reads spectrum and   
            headers.
        IndTable (actin2.IndTable) : Object containing the 
            indices table.
        ProcessSpec (actin2.ProcessSpec) : Object to process the 
            spectrum.
    """

    def __init__(self):
        self.ReadSpec = ReadSpec
        self.IndTable = IndTable
        self.CalcIndices = CalcIndices

    #@timeit
    def run(self, files, indices, obj_in=None, save_data=None, progress=False, verb=False, table_df=None, **spec_kw):
        """Run ACTIN for a list of files.

        Args:
            files (list) : List of fits files paths to be read.
            indices (list) : List of indices identification (as in 
                the indices table) to be calculated. If empty returns the fits headers only.
            rv_in (float, None) : If float, use this RV value to 
                shift spectrum to stellar rest frame, if None use RV value from fits file.
            obj_in (str, None) : If not None, the target name from 
                the fits file is overriden by this input.
            progress (bool) : If True show progress bar.
            table_df (None, pd.DataFrame) : Table with the indices parameters. If None use the built-in table.
            verb (bool) : Activate verbose option.
        
        Returns:
            pandas.DataFrame : Output pandas table.
        """

        printif("\nRUNNING ACTIN 2", verb)
        printif("---------------", verb)

        # Enforce 'files' type as list:
        if not isinstance(files, (list, np.ndarray, tuple)):
            files = [files]

        assert len(files) != 0, "'files' list is empty"

        files.sort()

        if progress:
            verb = False
            files = tqdm.tqdm(files)

        actin = [] 
        for file in files:
            time.sleep(0.01)
            read_spec = self.ReadSpec(file, verb=verb, **spec_kw)

            if not indices:
                return pd.DataFrame([read_spec.spec.headers])

            # read_spec.spec is the spectrograph class
            calc_ind = self.CalcIndices(read_spec.spec, table_df=table_df, indices=indices, interp=True, verb=verb)

            data = calc_ind.data

            data['actin_ver'] = '2.0.0'

            actin.append(data)


        df = pd.DataFrame(actin)

        # TODO: option to save rdb files

        if save_data:
            df.to_csv(save_data, index=False)
            printif(f"Data saved to: {save_data}", verb)

        return df

