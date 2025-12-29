function ptycho = a_loop_control_parameters()

% This paramfile will provide the parameters to 
% process the small simulated graphene example 4D STEM dataset
% which can be downloaded from 
% https://zenodo.org/record/3578284/files/m.mat?download=1

% This file is used to upload the initial parameters to perform the
% ptychography phase reconstructions. It creates a structured array named 'ptycho',
% which contains the intial variables for the functions to carry on
% calculating the different steps of the ptychography reconstruction. 
%
%
% Data identifier to know which is the data that has been analysed
ptycho.identifier = '';
% Acceleration voltage in [kV]
ptycho.voltage_kV  = 60;
% Wavelength in [Angtrom]
ptycho.wavelength = Wavelength(ptycho.voltage_kV); 
% Probe convergence angle in [Radian]
ptycho.ObjApt_angle = 20.0e-3; 
% Rotation angle in [Radian]
ptycho.rot_angle = 0;  
% Image pixel size in [Angstrom]. This corresponds to the sampling of the
% probe positions and it should be defined as [pixel_size_in_x pixel_size_in_y];
%ptycho.pix_size = [-(27.06-29.52)/25 -(21.3-25.56)/42]; % original
ptycho.pix_size = [0.1,0.1]; % zyd added

% Define which type of data is going to be loaded
% Flags are:
%       'PN' : for data acquired using PNdetector
%       'Sim': for data that was simulated using software
%	'mib': for data in mib format (Merlin)
%       'tif': for data that was acquired in tif
%       'hdf5': for data using HDF5 format 
%       'hdf5_mspy': for hdf5 data from mspy(0.3.4+)
ptycho.dataformat = '';
%
% If using simulation or tif, give the filename of the data
% Note that this can be a full path plus filename 
% or if you are running in the directory of the data, it can be just the filename 
%ptycho.data_filename = '../MSpy/MSpy-0.3.4_Graphene_rectangular_cell/test1_Graphene_rectangular_cell_CA25/dp_mat/dps.hdf5';
% ptycho.data_filename = '../mustem/dps_grc_padding_4dstem_01.hdf5';

% If using simulations, usually only one unit cell is simulated. Then, it
% can be repeated several times when the material is periodic. This flag
% determines if the simulation is going to be repeated. Value 1 (yes) and 0
% (no)
ptycho.repeat_unit_cell_simulation = 0;
% If unit cell is repeated, this value defines how many unit cells is
% repeated in x direction
ptycho.repeat_number_unit_cell_in_x = 5;
% If unit cell is repeated, this value defines how many unit cells is
% repeated in y direction
ptycho.repeat_number_unit_cell_in_y = 3;

% If one wishes to add poisson noise (e.g. to simulated data)
% set next line to 1 and specify dose in electrons per Angstrom squared
ptycho.add_noise = 0;
ptycho.dose = 1e5; % in e/A2

% Define if zero-padding is used. It is helpful to do this to avoid
% artifacts when doing the Fourier Transforms. Value 1 (yes) and 0
% (no). For experimental data it is usually recomended, for periodic
% simulated data it is not necessary.
ptycho.use_padding = 0;

% Flag to define if Dscan needs to be corrected. Values 1 (yes) or 0 (no)
% This function measures the centre of mass (COM) of the BF disc to quantify 
% any off-centreing of the BF disc due to D-Scan issue, which
% usually happens at low resolution without doing D-scan alignment tuning.
ptycho.IfCorrectDscan = 0;
% Flag to define is detector/ronchigrams were binned. Values 1 (binning) or
% 0 (no binning)
ptycho.detbinningFlag = 0;
% detbinningFlag = 0 will leave the detector as it is.
% detbinningFlag = 1 will bin the detector to the same pixel size. for example,
% if the detector is oiginally (132x264). detbinningFlag = 1
% bins the detector to (132x132)
% detbinningFlag = 4 will bin the detector to (66x66)
% detbinningFlag = 8 will bin the detector to (33x33)
%


% Flag to set if figures of each function are going to be plotted. Value 1
% (yes) and 0 (no)
ptycho.plot_figures = 0;
% Define value to use as threshold to select the bright field disk area in
% pacbed
ptycho.pacbedthreshold = 0.5; 
%
% For Wigner Distribution Deconvolution method, the probe aberrations can
% can be initialised here.
ptycho.aberr.C1   = 0  ; % Defocus      Unit: meters. (common value = 1nm)
ptycho.aberr.C12a = 0  ; % 2-fold Stig, a direction Unit: meters. (Common value = 5nm)
ptycho.aberr.C12b = 0  ; % 2-fold Stig, b direction Unit: meters. (Common value = 5nm)
ptycho.aberr.C21a = 0  ; % Axial Coma   Unit: meters. (Common value = 250nm)
ptycho.aberr.C21b = 0  ; % Axial Coma   Unit: meters. (Common value = 250nm)
ptycho.aberr.C23a = 0  ; % 3-fold Stig  Unit: meters. (Common value = 250nm)
ptycho.aberr.C23b = 0  ; % 3-fold Stig  Unit: meters. (Common value = 250nm)
ptycho.aberr.C30a = 0  ; % Cs           Unit: meters. (Common value = 50um)
ptycho.aberr.C32a = 0  ; % 3rd-order 2-fold Stig  Unit: meters. (Common value = 20um)
ptycho.aberr.C32b = 0  ; % 3rd-order 2-fold Stig  Unit: meters. (Common value = 20um)
ptycho.aberr.C34a = 0  ; % 4-fold Stig  Unit: meters. (Common value = 25um)
ptycho.aberr.C34b = 0  ; % 4-fold Stig  Unit: meters. (Common value = 25um)
ptycho.aberr.C50a = 0  ; % 5th ord. Cs  Unit: meters. (Common value = 50mm)
ptycho.aberr.C56a = 0  ; % 6-fold Stig  Unit: meters. (Common value = 50mm)
%
% epsilon ratio for WDD
ptycho.eps_ratio = 0.1;
ptycho.ChrDeconvFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varfunctions structure indicates which functions have been executed
ptycho.varfunctions.load_parameters = 0;
ptycho.varfunctions.load_data = 0;
ptycho.varfunctions.define_real_space = 0;
ptycho.varfunctions.detector_binning = 0;
ptycho.varfunctions.define_reciprocal_space = 0;
ptycho.varfunctions.calculate_detector_masks = 0;
ptycho.varfunctions.show_synthetic_IBF = 0;
ptycho.varfunctions.show_synthetic_BF = 0;
ptycho.varfunctions.show_synthetic_ABF = 0;
ptycho.varfunctions.show_synthetic_DF = 0;
ptycho.varfunctions.show_synthetic_DPC = 0;
ptycho.varfunctions.calculate_center_of_mass = 0;
ptycho.varfunctions.truncate_ronchigrams_within_collection_angle = 0;
ptycho.varfunctions.define_reciprocal_space_for_truncated_detector_plane = 0;
ptycho.varfunctions.FT_from_M_to_G = 0;
ptycho.varfunctions.single_side_band_reconstruction = 0;
ptycho.varfunctions.define_reciprocal_space_resampled_ronchigram = 0;
ptycho.varfunctions.define_probe_function = 0;
ptycho.varfunctions.FT_from_G_to_H = 0;
ptycho.varfunctions.wigner_distribution_deconvolution_reconstruction = 0;
ptycho.varfunctions.calculate_G_power_spectrum_wrt_Kf = 0;
% end of parameter file

end