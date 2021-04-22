"""
Test ACTIN 2
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import glob

import actin2

import glob

# TODO: Why are the indices different for HARPS s1d and e2ds files (CaII and Halpha?
# TODO: test the difference between applying interp or not in the indices
#? Ha16 is more appropriate for active M dwarfs - c.f. Gl551 - make tutorial to show this in readthedocs! 

#********************************
#* CONTROL PANEL:
test_spec = "HARPS"
test_file = "1d"

indices = ['I_CaII']
#indices = []



#********************************
#* Get files:

path = f'test/{test_spec}'

if test_spec == 'HARPS':
    if test_spec == "HARPS" and test_file == '1d':
        file_type = "_s1d_A"
    if test_spec == "HARPS" and test_file == '2d':
        file_type = "_e2ds_A"
            
    files = glob.glob(os.path.join(path, f"Gl551/*{file_type}.fits"))

if test_spec == 'HARPS-N':
    if test_spec == "HARPS-N" and test_file == '1d':
        file_type = "_s1d_A"
    if test_spec == "HARPS-N" and test_file == '2d':
        file_type = "_e2ds_A"

    files = glob.glob(os.path.join(path, f"*{file_type}.fits"))

if test_spec == 'ESPRESSO':
    if test_file == '1d':
        file_type = 'S1D_A'
    if test_file == '2d':
        file_type = 'S2D_A'

    files = glob.glob(os.path.join(path, f"*{file_type}.fits"))



files = np.sort(files)

#********************************

# if "HARPS-N" in test_spec:
#     path = 'test/HARPS-N'
#     if test_file == '1d':
#         files = glob.glob(os.path.join(path, "*_s1d_A.fits"))
#     if test_file == '2d':
#         files = glob.glob(os.path.join(path, "*_e2ds_A.fits"))

# if "ESPRESSSO" in test_mode:
#     path = "test/ESPRESSO"
#*
# path = "/Users/jgsilva/Astrophysics/2019_Ambre/Test_GP/HD41248"
# #* HARPS:
# #path = 'test/HARPS'
# #* HARPS-N:
# #path = "/Users/jgsilva/Astrophysics_blackdisk/2017_ActivityTool/fits_HARPSN_2d/HD80606"
# #path = 'test/HARPS-N'

# files = glob.glob(os.path.join(path, "*_s1d_A.fits"))
# #files = glob.glob(os.path.join(path, "*_e2ds_A.fits"))
# #* ESPRESSO:
# path = "/Users/jgsilva/Astrophysics/2019_Ambre/Test_GP"
# files = glob.glob(os.path.join(path, "HD41248_es/reduced/*_S1D_A.fits"))
# #files = glob.glob(os.path.join(path, "HD41248_es/reduced/*_S2D_A.fits"))
# path = "test/ESPRESSO"
# files = glob.glob(os.path.join(path, "*_S1D_A.fits"))
# #files = glob.glob(os.path.join(path, "*_S1D_FINAL_A.fits"))
# files = glob.glob(os.path.join(path, "*_S2D_A.fits"))
# #* SPIRou:
# path = 'test/SPIRou'
# files = glob.glob(os.path.join(path, "*s.fits"))

#********************************

#* Compare with ACTIN 1.3.6:
#from actin.actin import actin

# from decorators import timeit
# @timeit
# def run():
#     actin(files[:100], calc_index=[index])
# run()
# sys.exit()
# e2ds 100: very very slow due to the search of calib files!!!
#* s1d 100: 12.6284 sec


# output = os.path.join(sys.path[0], "ACTIN")
#* Compare with ACTIN 1.4:
# import actin2_old
# df = actin2_old.run_folder(files[0], index_id=[index], rv_in=None, obj=None, save_data=False, output_dir=None, bar=True, verb=True)
# # # #print(df)
# sys.exit()
#* e2ds 100: 11.6426 sec
#* s1d 100: 8.7248 sec

#* e2ds all: 30.2129 sec
#* s1d all: 21.3785 sec

#print(df[['rv', 'flg']])
#? ACTIN 1.4 is ~38% faster than ACTIN 1.3

#sys.exit()

#********************************
#* RUN

from actin2 import ACTIN

actin = ACTIN()

#* Use ACTIN to obtain the indices table:
#df = actin.IndTable().table



#* Use ACTIN to read just a spectrum:
#read_spec = actin.ReadSpec(files[0], verb=True)

# plot the spectrum
#read_spec.plot(order=56, show=True, c='b')
#read_spec.plot_ccf(show=True)
#read_spec.plot_bisector(show=True)

# spec = read_spec.spec
# headers = spec.headers
# spectrum = spec.spectrum

#print(spec.ccf_bisector.keys())
#print(pd.DataFrame([headers]))

#sys.exit()


#plt.errorbar(spectrum['wave'], spectrum['flux'], spectrum['flux_err'], fmt='k.')
#plt.show()

#* Calculate indices:
#calc_ind = actin.CalcIndices(spec, indices, full_output=False)
#data = calc_ind.data

#print(data.keys())

# print()
# if indices:
#     print(data[indices[0]], data[indices[0] + "_err"])
# else:
#     print(data)

#* To test the spec_class_in
from spectrographs.HARPN import HARPN
spec_in = HARPN


#sys.exit()

df = actin.run(files[0], indices, verb=False, progress=False, obj_in=None, get_bis=True, get_ccf=True, spec_class_in=None)
print(df.keys())
print(df)

# if indices:
#     print(df[["rv", "rv_err", indices[0], indices[0] + "_err", 'spec_flg']])
# else:
#     print(df)

#print(df[indices[0]].mean())

#print(df[['obj', 'instr', 'ftype', 'date_obs', 'bjd', index, index + '_err', index + '_flg', 'spec_flg', 'ron', 'rv', 'fwhm', 'bis', 'berv', 'seeing_start']])
#print(df[['obj', 'instr', 'date_obs', 'bjd', index, index + '_err', index + '_flg', 'spec_flg', 'ron', 'ccf_noise', 'rv_err', 'airmass', 'seeing']])

#* e2ds 100: 15.1674 sec
#* s1d 100: 14.1615 sec

#* e2ds all: 33.4152 sec
#* s1d all: 29.5886 sec


#!
# EXAMPLE OF INHERITING FROM A CLASS:
# import actin2

# class preprocess(actin2.ReadSpec):

#     def __init__(self, file, instr, rv_in=None, sep_fiber=True, ccf=True):
#         self.spec = actin2.ReadSpec(file, instr).spec


#     def plot_spec(self):
#         plt.plot(self.spec.spectrum['wave'][7], self.spec.spectrum['flux'][7])
#         plt.show()



# spec = preprocess(files[0], instr='HARPS')
# spec.plot_spec()
# sys.exit()
#!