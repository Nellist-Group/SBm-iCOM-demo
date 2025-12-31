# Demostration of SBm-iCOM on Matlab

This is the code demostration for article ___Defocus correction and noise reduction using a hybrid ptychography and Centre-of-Mass algorithm___ (DOI: 10.1111/jmi.70010). The dataset for the demo is available at: https://zenodo.org/records/15011384 . The demo was tested on MATLAB 2020b.

updated on 31st Dec 2025.

## To run the demo:

1. Download demonstration dataset from https://zenodo.org/records/15011384.
2. Change `dataPath` (line 2) in the scirpt `Simu4dstem_script_for_exp_data_df0.m` to the path of dataset file `binned_diff_20230819_163033_df0.hdf5`. For example:

  ```
  % set experimental parameters
  dataPath = './the path to binned_diff_20230819_163033_df0.hdf5'; % datapath for data file
  ```
3. Change `savePath` (line 5) in the script to where the results need to be saved. For example:

  ```
  savePath = './some folder to save results/';
  ```

4. save the script and run it.
5. change `dataPath` for scirpt `Simu4dstem_script_for_exp_data_df10.m` to dataset file `binned_diff_20230819_163214_df10.hdf5`, and script `Simu4dstem_script_for_exp_data_df20.m` to dataset file `binned_diff_20230819_163728_df20.hdf5` . Edit `savePath` if needed, then save them and run them.
6. check results in savePath.
   | result file name | description |
   | ---- | ---- |
   | result_img/img_COMx.png | COMx without rotation correction |
   | result_img/img_COMy.png | COMy without rotation correction |
   | result_img/img_COMxf.png | COMx without rotation correction, high-pass filtered |
   | result_img/img_COyxf.png | COMy without rotation correction, high-pass filtered |
   | result_img/img_COMxr.png | COMx with rotation correction |
   | result_img/img_COMyr.png | COMy with rotation correction |
   | result_img/img_iCOM.png | iCOM without high-pass filtering |
   | result_img/img_iCOMf.png | iCOM, high-pass filtered |
   | result_img/img_dCOM.png | dCOM |
   | result_img/img_SBmCOMx.png | SBmCOMx without rotation correction |
   | result_img/img_SBmCOMy.png | SBmCOMy without rotation correction |
   | result_img/img_SBmCOMxr.png | SBmCOMx with rotation correction |
   | result_img/img_SBmCOMyr.png | SBmCOMy with rotation correction |
   | result_img/img_SBmiCOM.png | SBmiCOM without high-pass filtering |
   | result_img/img_SBmiCOMf.png | SBmiCOM, high-pass filtered |
   | result_img/img_SBmdCOM.png | SBmdCOM |
   | result_img/img_SSBAC_amp_l.png | SSBAC amplitude, left side-band |
   | result_img/img_SSBAC_amp_r.png | SSBAC amplitude, right side-band |
   | result_img/img_SSBAC_phs_l.png | SSBAC phase, left side-band |
   | result_img/img_SSBAC_phs_r.png | SSBAC phase, right side-band |
   | result_img/img_STEM_25_50.png | integrated ADF, collection angle 25 to 50 mrad |
   | result_img/img_STEM_60_200.png | integrated HAADF, collection angle 60 to 200 mrad |
   | result_mat/COM.mat | COM results in mat file |
   | result_mat/SBmCOM.mat | SBmCOM results in mat file |
   | result_mat/SSBAC.mat | SSBAC parameters and results in mat file |
   | result_mat/STEM.mat | STEM results (ADF and HAADF) in mat file |

## Reference

COM/DPC:

getDPC (https://github.com/hachteja/GetDPC/tree/master)

ptychoSTEM package:

T. J. Pennycook, A. R. Lupini, H. Yang, M. F. Murfitt, L. Jones, P. D. Nellist, Ultramicroscopy 151 (2015) 160 – 167.

H. Yang, T.J. Pennycook, P.D. Nellist, Ultramicroscopy, 151 (2015) 232-239.

H. Yang, R.R. Rutte, L. Jones, M. Simson, R. Sagawa, H. Ryll, M. Huth, T.J. Pennycook, H. Soltau, Y. Kondo, B. Davis. Nature Communications 7 (2016) 12532.
