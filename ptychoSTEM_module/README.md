## ptychoSTEM

ptychoSTEM is an open source set of MATLAB scripts that can be used to perform non-iterative ptychographic phase reconstructions with 4D STEM data. Currently the code implements the Single Side Band (SSB) method [1,2] and the Wigner Distribution Deconvolution (WDD) method [3].

The code is based on a set of functions that use a structured array called ptycho. This array contains all the variables needed for each function and each function will modify or generate new variables inside ptycho. Many of the functions have comments that explain their purpose and requirements. 

The files SSB_main.m and WDD_main.m contain the sequence of functions needed to perform either SSB or WDD ptychography, with comments to explan the process. A parameter file is used to contain all the parameters used for processing and an example is shown in the file example_parameter_file.m, also with some comments to help explain the parameters.

A small (~64 MB) simulated graphene example 4D STEM dataset can be downloaded here: https://zenodo.org/record/3578284/files/m.mat?download=1. 
To process this example dataset you can use the simulted_graphene_example_paramfile.m file contained in the ptychoSTEM repository 
as the parameter file. It provides parameters that should work for the example dataset. 

## Running the code

A suggested workflow is to first load the ptychoSTEM repository into MATLAB's path, then change to the directory containing the data to process. If the parameter file is named paramfile.m and placed in the directory containing the data it will be loaded automatically when running the code. If you use a separate directory for each dataset the paramfile then can also serve as a record of what parameters you used to process the data. This strategy also generally keeps code separate from parameters. 

For example if you download the above mentioned small graphene dataset to a directory/folder then copy from the repository simulted_graphene_example_paramfile.m to paramfile.m inside the same directory as the data, then the code will automatically load the correct paramters and process the dataset (assuming you are in the same directory as the data in Matlab and the code is in Matlab's path). 

Running the functions in SSB_main.m or WDD_main.m sequentially will go through the steps needed to process the data for SSB or (SSB and) WDD. If all the parameters are specified correctly in the parameter file, they can be run one after the other immediately either by selecting all the functions you want to run and performing a run selection in MATLAB or even just running the entire script. 

The file condensed_WDD_Main.m provides a similar sequence of funtions as WDD_main.m (including the functions for SSB) but in condensed form without many comments, which can be more convienent for actually running the code, especially once you have some idea of how things work.

#### Preparing the data
Running the following sequence of functions is fairly minimal preperation procedure
```matlab
[ptycho] = load_parameters();
[ptycho] = load_data(ptycho);
[ptycho] = define_real_space(ptycho);
[ptycho] = define_reciprocal_space(ptycho);
[ptycho] = truncate_ronchigrams_within_collection_angle(ptycho);
[ptycho] = define_reciprocal_space_for_truncated_detector_plane(ptycho);
[ptycho] = FT_from_M_to_G(ptycho);
```
this gets us to the point of having performed the Fourier transform with respect to probe position of the data, i.e. up to calculating G(K,Q) in reference 1. 
(Assuming we had some workable parameters set in a paramfile.m in the current directory.)

#### Calibrating the trotter selection

However **it is important that the calibration of the locations of the double-disk overlaps (AKA "the trotters" as they are colloquially termed) is correct**, and once we have reached the point of calculating G(K,Q) it is possible to check this by observing how the code is chopping out the trotters. An example of chopping out the double disk overlap regions for a particular spatial frequency (Q) can be seen by comparing Fig 3 with Fig 4 of reference 1. A tool that facilitates this is 
```matlab
    show_trotters(ptycho,5)
```
which will automatically choose the 5 frequencies with the strongest amplitude in the data (usually you will always have ptycho as the first variable, but you can change 5 to whatever number of frequencies you like), and display both the phase and amplitude of each of those frequencies over the pixels of the detector and how the double disk overlaps are being cut out with the current calibration. You can then adjust the rotation and other calibrations such that the double disk overlaps are best cut out, running show_trotters again after each adjustment to check how the fit improves (or not). 


The parameters that are important to the calibration of the trotters are the <br/>
* rotation angle (ptycho.rot_angle) <br/>
* probe convergence angle (ptycho.ObjApt_angle)<br/>
* acceleration voltage (ptycho.voltage_kV) which is converted into wavelength (ptycho.wavelength)<br/>
* probe position pixel size or scan step size (ptycho.pix_size)<br/>


Typically one sets these in the parameter file, and can be recorded during experiments with the common exception of the rotation angle. The rotation angle applies a uniform rotation of the calibration around the optical axis. This is usually not needed for simulations, but changes in magnetic lens settings in a microscope can cause such a rotation. Setting the rotation to a sensible value is critical to correctly selecting the disk overlaps. For weak objects it is often relatively straightforward to see where the double disk overlaps occur in probe reciprocal space. 

Note that if you need to adjust ptycho.pix_size you will want to rerun 
```matlab
[ptycho] = define_real_space(ptycho);
```

Once the calibration is correct you can proceed with the ptychographic processing. 

For SSB, one then just needs to run
```matlab
[ptycho] = single_side_band_reconstruction(ptycho);
```
and if automatic displaying is not enabled run
```matlab
show_SSB_results(ptycho)
```
to display the results in a figure.

For WDD, you will want to additionally run
```matlab
[ptycho] = define_reciprocal_space_resampled_ronchigram(ptycho);
[ptycho] = define_probe_function(ptycho);
[ptycho] = FT_from_G_to_H(ptycho);
[ptycho] = wigner_distribution_deconvolution_reconstruction(ptycho);
```
and similarly if automatic displaying is not enabled run
```matlab
show_WDD_results(ptycho)
```

#### Conventional modes from 4D datasets
Conventional STEM imaging can also be synthesized by integrating the diffraction patterns over a virtual detector geometry, including:
Bright Field (BF), Incoherent Bright Field (IBF), Annular Bright Field (ABF), Annular Dark Field (ADF) and the calculation of the Center of Mass (COM) of the bright field disk. These can be useful when debugging issues such as data stored in an unexpected order. To generate these signals you will want to run
```matlab
[ptycho] = calculate_detector_masks(ptycho);
```

after getting to the point of defining reciiprical space in the code ie after
```matlab
[ptycho] = define_reciprocal_space(ptycho);
```
in the previous "preparation" sequence of functions then for example to show ADF and ABF images run 

```matlab
[ptycho] = show_synthetic_DF(ptycho);
[ptycho] = show_synthetic_ABF(ptycho);
```


In case any publications result from using this software, please cite the following publications:<br/>
    - for the SSB method: <br/>
        *[1] [T. J. Pennycook, A. R. Lupini, H. Yang, M. F. Murfitt, L. Jones, P. D. Nellist, Ultramicroscopy 151 (2015) 160 – 167.](https://doi.org/10.1016/j.ultramic.2014.09.013)*<br/>
        *[2] [H. Yang, T.J. Pennycook, P.D. Nellist, Ultramicroscopy, 151 (2015) 232-239.](https://doi.org/10.1016/j.ultramic.2014.10.013)*<br/>
    - for the WDD method:<br/>
        *[3] [H. Yang, R.R. Rutte, L. Jones, M. Simson, R. Sagawa, H. Ryll, M. Huth, T.J. Pennycook, H. Soltau, Y. Kondo, B. Davis. Nature Communications 7 (2016) 12532.](https://doi.org/10.1038/ncomms12532)*<br/>

Further suggested reading:<br/>
*[4] [J.M. Rodenburg, et al. Ultramicroscopy 48 (1993) 304-314](https://doi.org/10.1016/0304-3991(93)90105-7)*.<br/>
*[5] [J.M. Rodenburg and R.H. Bates, Phil. Trans. R. Soc. Lond. A (1992) 339, 521-553](https://doi.org/10.1098/rsta.1992.0050)*<br/>
*[6] [P.D. Nellist, B.C. McCallum, J.M. Rodenburg, Nature, 374 (1995) 630-632.](https://doi.org/10.1038/374630a0)*









