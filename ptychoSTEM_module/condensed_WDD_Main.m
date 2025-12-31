[ptycho] = load_parameters();
[ptycho] = load_data(ptycho);
if ptycho.add_noise
    ptycho.m = add_noise_to_CBEDs(ptycho.m, ptycho.dose, ptycho.pix_size);
end
[ptycho] = define_real_space(ptycho);
if ptycho.detbinningFlag
    [ptycho] = detector_binning(ptycho) ;
end
[ptycho] = define_reciprocal_space(ptycho);
[ptycho] = calculate_detector_masks(ptycho);
[ptycho] = show_synthetic_DF(ptycho);
[ptycho] = show_synthetic_ABF(ptycho);


[ptycho] = truncate_ronchigrams_within_collection_angle(ptycho);

if 0 %use if the BF disk cropped already use this instead of previous line
    ptycho.det_range_x=1:size(ptycho.m,1);
    ptycho.det_range_y=1:size(ptycho.m,2);
    ptycho.m_wp=ptycho.m;
end
    
[ptycho] = define_reciprocal_space_for_truncated_detector_plane(ptycho);
[ptycho] = FT_from_M_to_G(ptycho);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    [ptycho] = FT_from_M_to_G(ptycho);
%show_trotters(ptycho,1)
%run to here for testing trotters

[ptycho] = single_side_band_reconstruction(ptycho);

[ptycho] = define_reciprocal_space_resampled_ronchigram(ptycho);
[ptycho] = define_probe_function(ptycho);
[ptycho] = FT_from_G_to_H(ptycho);
[ptycho] = wigner_distribution_deconvolution_reconstruction(ptycho);
%local_tif_saver

%save_output_tif
