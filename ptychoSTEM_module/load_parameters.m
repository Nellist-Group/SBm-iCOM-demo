function [ptycho]=load_parameters()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To load parameters. This function opens and runs a parameter m-file, which
% contains the definition of the starting variables included in the
% structured array named ptycho.
% It adds to ptycho:
% ptycho.varfunctions.load_parameters: Flag to indicate the function has
% been run
% ptycho.varfunctions.load_parameters_pathname: string containing the path
% where the parameter file was loaded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenameParameters='paramfile.m'; %Defaults
pathname = pwd;%full path to current directory
defaultfile=strcat(pathname,filesep,filenameParameters);

if exist(defaultfile, 'file')
  % Default Paramter File exists.
  disp('Using default parameter file')
  cd(pathname);
  run(filenameParameters);

else
  % Default File does not exist. Ask user.
 [filenameParameters, pathname] = uigetfile('*.m', 'Select user input parameters');
  cd(pathname);
  run(filenameParameters);
end

ptycho.varfunctions.load_parameters = 1;
ptycho.varfunctions.load_parameters_pathname = pathname;


end
