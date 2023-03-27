function get_new_Cmem(CELL_LEN, CELL_WID, CELL_DEP, population_size)
    Cmem_array  = [];
    for tub_no = 1:population_size
        Cmem = get_single_Cmem(CELL_LEN, CELL_WID, CELL_DEP, tub_no);
        Cmem_array = [Cmem_array, Cmem];
    end

    save('get_new_Cmem.mat');
end

function Cmem = get_single_Cmem(nx, ny, nz, tub_no)
        filename = ['./TRIAL/tub_input_ver2_' num2str(tub_no) '.txt'];
        data        = load(filename);
        num_tubule = sum(data);
        %%
        specific_cmem_surface    = 1e-2;         % [pF/um^2]
        specific_cmem_tubule     = 0.56e-2;      % [pF/um^2]

        cru_length       = 1.84;
        cru_width        = 0.9;
        pi               = 3.1415926;
        diameter_tubule  = 0.3;
        at_tt_ratio      = 6.0;
        length   = cru_length * nx;
        width    = cru_width * ny;
        depth    = cru_width * nz;

        area_surf        = 2.0 * (length * width + length * depth + width * depth);
        area_single_tub  = pi * diameter_tubule * ( at_tt_ratio/(at_tt_ratio + 2.0) * cru_length + 2.0/(at_tt_ratio + 2.0) * cru_width );
        num_surf_cru     = (nx * ny ) * 2 + ( (nz - 2) * nx ) * 2 + ( (ny - 2) * (nz - 2) )* 2;  % 3130 (one layer) rather than 5644 (two layers)

        % the most densely tubulated cell has 8998 coupled CRUs - Assumption Cmem = 110 [pF]
        % e.g., the de-tubulated cell has 3130 surface CRUs and (5644-3130) second-layer inner CRUs
        % num_inner_cru_dense_tub  = num_dense_tubulated_CRU * 1.0 - num_surf_cru;
        % area_sum_dense_tub       = num_inner_cru_dense_tub * area_single_tub + area_surf;
        num_inner_cru    = num_tubule - num_surf_cru;
        Cmem    = (area_surf * specific_cmem_surface + num_inner_cru * area_single_tub * specific_cmem_tubule) * 1e-12; % [F]
end