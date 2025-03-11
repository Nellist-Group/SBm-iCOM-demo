function [atoms,lx,ly,lz] = Simu4dstem_preprocess_atoms(jsonData)
%MULTEM_4DSTEM_PREPROCESS_ATOMS Summary of this function goes here
%   Detailed explanation goes here

% check para
if jsonData.nx ~= jsonData.ny
    error('Simu4dstem_preprocess_atoms: jsonData.nx ~= jsonData.ny');
end

load(jsonData.model_path, jsonData.model_variable)
eval(['atoms = ',jsonData.model_variable,';']);

atoms(:, 4) = atoms(:, 4) - min(atoms(:, 4));

% Temporary solution for thickness bug
for ii = size(atoms,1):-1:1
    if atoms(ii,4) > jsonData.thick_max
        atoms(ii,:) = [];
    end
end

box = num2cell(max(atoms(:, 2:4)));
[lx_atom, ly_atom, lz_atom] = deal(box{:});

extra_pixel_num = jsonData.extra_pxl_num_on_atom_pot; % extra pixels on the edge of calcualted atom pot.
lambda = 1.23984244/sqrt(jsonData.voltage*(2*510.99906+jsonData.voltage))*1e1;% A, calculate wave length
dxy = lambda / (jsonData.dk * 1e-3) / jsonData.nx; % in A
scan_size_A = [jsonData.step_size*(jsonData.array_size_x-1),jsonData.step_size*(jsonData.array_size_y-1)]; % in A, no probe mat width
probe_mat_width_A = [jsonData.nx,jsonData.ny] .* dxy; % in A
atom_pot_size_A = scan_size_A + probe_mat_width_A + 2.*extra_pixel_num.*dxy;
atom_pot_size_pxl = round(atom_pot_size_A ./ dxy);
lx = atom_pot_size_A(1);
ly = atom_pot_size_A(2);

atoms(:,2) = atoms(:,2) + (extra_pixel_num.*dxy + (jsonData.nx/2).*dxy - jsonData.first_xy(1));
atoms(:,3) = atoms(:,3) + (extra_pixel_num.*dxy + (jsonData.ny/2).*dxy - jsonData.first_xy(2));

for ii = size(atoms,1):-1:1
    if atoms(ii,2) < 0 || atoms(ii,3) < 0
        atoms(ii,:) = [];
    end
end

% calculate lx,ly by dk setting
% lx = lambda/(jsonData.dk*1e-3);
% ly = lx;
lz = lz_atom;

end

