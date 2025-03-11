function ptycho = a_loop_control_read4ddataset(Im_Path, ptycho, var_name)
%A_LOOP_CONTROL_READ4DDATASET Summary of this function goes here
%   input: varname --> the name of the variable in the mat file

EXT = Im_Path(end-3:end);

Done = 0;
switch EXT
    case '.mat'
        load(Im_Path,var_name);
        eval(['m = ',var_name,';']);
        %%%%%%%% add for multem simulated datasets %%%%%%%%
        clear(var_name);
        m = permute(m,[1,2,4,3]);
        %%%%%%%% end for multem simulated datasets %%%%%%%%
        if exist('m','var') && size(m,4) > 1
            ptycho.m = m;
            clear m
            Done = 1;
        else
            ptycho = rmfield(ptycho,'m');
            error('".mat" file does not contain a 4D "m" matrix.');
        end
        
    case 'hdf5'
        hinfo = h5info(Im_Path);
        if isfield(hinfo.Groups.Groups,'Name')
            ptycho.m = double(h5read(Im_Path,horzcat(hinfo.Groups.Groups.Name,'/data')));
            if size(ptycho.m,4) > 1
                Done = 1;
            else
                ptycho = rmfield(ptycho,'m');
                error( 'Could not find 4D dataset.');
            end
        end
        
        %%%%%%%%%%% Zhiyuan Ding added %%%%%%%%%%%
        % load hdf5 data from muSTEM2mat
    case 'muh5'
        ptycho.m = h5read(Im_Path,'/dps');
        Done = 1;
        
    case 'frh5'
        % hdf5 file from frsm6_convertor
        ptycho.m = h5read(Im_Path,'/dps');
        m_dim = size(ptycho.m);
        if length(m_dim) == 3
            total_dp_num_sqrt = sqrt(m_dim(3));
            if total_dp_num_sqrt == fix(total_dp_num_sqrt)
                ptycho.m = reshape(ptycho.m,[m_dim(1),m_dim(2),total_dp_num_sqrt,total_dp_num_sqrt]);
            end
        end
        ptycho.m = single(ptycho.m);
        Done = 1;
        %%%%%%%%%% Zhiyuan Ding added end %%%%%%%%%%%
end
if Done == 1
    if (~isa(ptycho.m,'double')) && (~isa(ptycho.m,'single'))
        ptycho.m = single(ptycho.m);
    end
    ptycho.probe_range = [0 0 size(ptycho.m,4) size(ptycho.m,3)];
    ptycho.varfunctions.load_data = 1;
end

end
