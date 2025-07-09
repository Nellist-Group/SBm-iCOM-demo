# Demostration of SBm-iCOM on Matlab

This is the code demostration for article ___Defocus correction and noise reduction using a hybrid ptychography and Centre-of-Mass algorithm___ (DOI: 10.1111/jmi.70010). The dataset for the demo is available at: https://zenodo.org/records/15011384 . The demo was tested on MATLAB 2020b.

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

## Reference

COM/DPC:

getDPC (https://github.com/hachteja/GetDPC/tree/master)

ptychoSTEM package:

T. J. Pennycook, A. R. Lupini, H. Yang, M. F. Murfitt, L. Jones, P. D. Nellist, Ultramicroscopy 151 (2015) 160 – 167.

H. Yang, T.J. Pennycook, P.D. Nellist, Ultramicroscopy, 151 (2015) 232-239.

H. Yang, R.R. Rutte, L. Jones, M. Simson, R. Sagawa, H. Ryll, M. Huth, T.J. Pennycook, H. Soltau, Y. Kondo, B. Davis. Nature Communications 7 (2016) 12532.
