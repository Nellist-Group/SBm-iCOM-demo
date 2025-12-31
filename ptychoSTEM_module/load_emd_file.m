function s = load_emd_file(f)
% LOAD_EMD_FILE Load a EMD file
%   s = load_emd_file('a_file.emd')
%
%   Returns a struct with all the signals contained
%   in the file, which is all the groups having an
%   attribute emd_group_type: 1, and a dataset called "data".
%
%   Currently it loads no metadata, and is not heavily tested.
%   So there is a high probability for it not working for most
%   EMD-files.

    fid = H5F.open(f);
    hdf5info = h5info(f);
    emd_group_list = findemdgroups(hdf5info);
    emd_group_list = emd_group_list{1};
    s = struct();
    for i=1:length(emd_group_list)
        emd_group = emd_group_list{i};
        signal_struct = struct('data',h5read(f,[emd_group '/data']));
        data_dim = size(signal_struct.data);
        gid = H5G.open(fid,emd_group);
        for j=1:length(data_dim)
            dim_tag = ['dim' num2str(j)];
            if H5L.exists(gid,dim_tag,'H5P_DEFAULT')
                signal_struct.(dim_tag) = h5read(f,[emd_group '/' dim_tag]);
            else
                signal_struct.(dim_tag) = 0:data_dim(j)-1;
            end
        end
        signal_tag = ['signal' num2str(i)];
        s.(signal_tag) = signal_struct;
    end
end

function emd_group_list = findemdgroups(hdf5struct)
    emd_group_list = {};
    emd_group_list = find_emd_attribute_in_groups(hdf5struct, emd_group_list);
end

function string_list = find_emd_attribute_in_groups(hdf5struct, string_list)
    group_list = getfield(hdf5struct,'Groups');
    for i = 1:length(group_list)
        group = group_list(i);
        if not(isempty(group.Attributes))
            bool = isemdgroup(group.Attributes);
            if bool
                string_list{end+1} = group.Name;
            end
        end
        if not(isempty(group))
            temp_list = findemdgroups(group);
            if not(isempty(temp_list))
                string_list{end+1} = temp_list;
            end
        end
    end
end

function bool = isemdgroup(attributes)
    for i = 1:length(attributes)
        bool = strcmp(attributes(i).Name,'emd_group_type');
        if bool
            break
        end
    end
end
