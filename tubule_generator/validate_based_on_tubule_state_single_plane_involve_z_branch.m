% To validate the tubule properties based on tubule structures
function validate_based_on_tubule_state_single_plane_involve_z_branch(TUBULE_POPULATION_NUM, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN)

    % 0: z-stacks average
    % 1: central plane
    % 2: bottom plane
    flag_palne  = 1;
    
    validate_tubule_at_branch       = zeros(TUBULE_POPULATION_NUM, CELL_LEN * CELL_WID * CELL_DEP);
    validate_tubule_tt_branch       = zeros(TUBULE_POPULATION_NUM, CELL_LEN * CELL_WID * CELL_DEP);
    validate_tubule_z_tt_branch     = zeros(TUBULE_POPULATION_NUM, CELL_LEN * CELL_WID * CELL_DEP);

    validate_tubule_dist_means      = zeros(TUBULE_POPULATION_NUM, 1);
    validate_tubule_dist_stds       = zeros(TUBULE_POPULATION_NUM, 1);
    validate_tubule_dist_medians    = zeros(TUBULE_POPULATION_NUM, 1);
    
    validate_cru_tubule_in_plane    = zeros(TUBULE_POPULATION_NUM, 3);
    for id_file = 1 : TUBULE_POPULATION_NUM
        [validate_tubule_at_branch(id_file, :), validate_tubule_tt_branch(id_file, :), validate_cru_tubule_in_plane(id_file, :), validate_tubule_z_tt_branch(id_file, :)]   = validate_tubule_collect_branch_length_by_branch(flag_palne, id_file, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER);
        [validate_tubule_dist_means(id_file, :), validate_tubule_dist_stds(id_file, :), validate_tubule_dist_medians(id_file, :)]                                           = validate_tubule_collect_inter_distance(id_file, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER);
    end
    validate_tubule_at_length           = sum(validate_tubule_at_branch')';
    validate_tubule_tt_length           = sum(validate_tubule_tt_branch')';
    validate_tubule_at_tt               = validate_tubule_at_length ./ validate_tubule_tt_length;
    validate_tubule_at_tt(isnan(validate_tubule_at_tt))     = 0;
    validate_tubule_at_tt_mean          = mean(validate_tubule_at_tt(validate_tubule_at_tt ~= 0));
    validate_tubule_at_tt_std           = std(validate_tubule_at_tt(validate_tubule_at_tt ~= 0));
    validate_tubule_branch_length_mean  = mean([validate_tubule_at_branch(validate_tubule_at_branch>0)' validate_tubule_tt_branch(validate_tubule_tt_branch>0)']);
    validate_tubule_branch_length_std   = std([validate_tubule_at_branch(validate_tubule_at_branch>0)' validate_tubule_tt_branch(validate_tubule_tt_branch>0)']);
    
    validate_tubule_branch_length_involve_z_dot_mean    = mean([validate_tubule_at_branch(validate_tubule_at_branch>0)' validate_tubule_tt_branch(validate_tubule_tt_branch>0)' validate_tubule_z_tt_branch(validate_tubule_z_tt_branch>0)']);
    validate_tubule_branch_length_involve_z_dot_std     = std([validate_tubule_at_branch(validate_tubule_at_branch>0)' validate_tubule_tt_branch(validate_tubule_tt_branch>0)' validate_tubule_z_tt_branch(validate_tubule_z_tt_branch>0)']);
    
    [validate_tubule_t_index_inner, validate_tubule_t_index_whole, validate_tubule_density_inner, validate_tubule_density_whole]    = validate_tubule_calculate_density_t_index(flag_palne, validate_tubule_at_length, validate_tubule_tt_length, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN);
    
    validate_tubule_tt_length_involve_z_dot     = sum(validate_tubule_tt_branch')' + sum(validate_tubule_z_tt_branch')';
    [validate_tubule_t_index_inner_involve_z_dot, validate_tubule_t_index_whole_involve_z_dot, validate_tubule_density_inner_involve_z_dot, validate_tubule_density_whole_involve_z_dot]    = validate_tubule_calculate_density_t_index(flag_palne, validate_tubule_at_length, validate_tubule_tt_length_involve_z_dot, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN);

    validate_tubule_at_tt_involve_z_dot = validate_tubule_at_length ./ validate_tubule_tt_length_involve_z_dot;
    validate_tubule_at_tt_involve_z_dot(isnan(validate_tubule_at_tt_involve_z_dot))     = 0;
    validate_tubule_at_tt_mean_involve_z_dot          = mean(validate_tubule_at_tt_involve_z_dot(validate_tubule_at_tt_involve_z_dot ~= 0));
    validate_tubule_at_tt_std_involve_z_dot           = std(validate_tubule_at_tt_involve_z_dot(validate_tubule_at_tt_involve_z_dot ~= 0));
    
    validate_tubule_density_inner_mean  = mean(validate_tubule_density_inner);
    validate_tubule_density_inner_std   = std(validate_tubule_density_inner);
    validate_tubule_density_whole_mean  = mean(validate_tubule_density_whole);
    validate_tubule_density_whole_std   = std(validate_tubule_density_whole);
    
    validate_tubule_t_index_inner_mean  = mean(validate_tubule_t_index_inner);
    validate_tubule_t_index_inner_std   = std(validate_tubule_t_index_inner);
    validate_tubule_t_index_whole_mean  = mean(validate_tubule_t_index_whole);
    validate_tubule_t_index_whole_std   = std(validate_tubule_t_index_whole);
    
    validate_tubule_dist_mean    = mean(validate_tubule_dist_means);
    validate_tubule_dist_std     = std(validate_tubule_dist_stds);
    
    validate_tubule_density_inner_involve_z_dot_mean  = mean(validate_tubule_density_inner_involve_z_dot);
    validate_tubule_density_inner_involve_z_dot_std   = std(validate_tubule_density_inner_involve_z_dot);
    validate_tubule_t_index_inner_involve_z_dot_mean  = mean(validate_tubule_t_index_inner_involve_z_dot);
    validate_tubule_t_index_inner_involve_z_dot_std   = std(validate_tubule_t_index_inner_involve_z_dot);    
    
    save(['validation_based_on_tubule_structure_single_plane_involve_z_branch.mat']) 
end

function [validate_tubule_at_branch, validate_tubule_tt_branch, validate_cru_tubule_in_plane, validate_tubule_z_tt_branch]     = validate_tubule_collect_branch_length_by_branch(flag_palne, id_file, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER)
    load(['TRIAL/simulation_log_' num2str(id_file) '.mat']);
    validate_tubule_at_branch  = zeros(1, CELL_LEN * CELL_WID * CELL_DEP);
    validate_tubule_tt_branch  = zeros(1, CELL_LEN * CELL_WID * CELL_DEP);
    branch_at_no_tmp    = 1;
    branch_tt_no_tmp    = 1;
    
    cru_x_branch        = 0;
    cru_y_branch        = 0;
    cru_z_branch        = 0;
    
    if (flag_palne == 0)
        for z_tmp = (1 + SURFACE_LAYER) : (CELL_DEP - SURFACE_LAYER)
    %     z_tmp = int64(CELL_DEP / 2);
    %     z_tmp = (1 + SURFACE_LAYER);
            % To filter all branches (>=2CRU in a certain xy-plate)
            tubule_cru_coordiante = tubule_cru_coordinate((tubule_cru_coordinate(:, 3) == z_tmp), :);
            for cru_id = 1 : (length(tubule_cru_coordiante(:, 1)) - 1)
                if ( (tubule_cru_coordiante(cru_id, 1)==tubule_cru_coordiante(cru_id+1, 1)) + (tubule_cru_coordiante(cru_id, 2)==tubule_cru_coordiante(cru_id+1, 2)) == 1)
                    validate_tubule_at_branch(1, branch_at_no_tmp)  = validate_tubule_at_branch(1, branch_at_no_tmp) + (tubule_cru_coordiante(cru_id, 1)~=tubule_cru_coordiante(cru_id+1, 1)) * CRU_LEN;
                    validate_tubule_tt_branch(1, branch_tt_no_tmp)  = validate_tubule_tt_branch(1, branch_tt_no_tmp) + (tubule_cru_coordiante(cru_id, 2)~=tubule_cru_coordiante(cru_id+1, 2)) * CRU_WID;
                else
                    branch_at_no_tmp    = branch_at_no_tmp + (validate_tubule_at_branch(1, branch_at_no_tmp) ~= 0);
                    branch_tt_no_tmp    = branch_tt_no_tmp + (validate_tubule_tt_branch(1, branch_tt_no_tmp) ~= 0);
                end
            end        
            branch_at_no_tmp    = branch_at_no_tmp + 1;
            branch_tt_no_tmp    = branch_tt_no_tmp + 1;

            for cru_id = 1 : (length(tubule_cru_coordiante(:, 1)) - 1)
                if(tubule_cru_coordiante(cru_id, 1)>=(1 + SURFACE_LAYER) && tubule_cru_coordiante(cru_id, 1)<=(CELL_LEN - SURFACE_LAYER))
                    if(tubule_cru_coordiante(cru_id, 2)>=(1 + SURFACE_LAYER) && tubule_cru_coordiante(cru_id, 2)<=(CELL_WID - SURFACE_LAYER))
                        cru_z_branch        = cru_z_branch + 1;
                    end
                end
            end
        end
    else
        if (flag_palne == 1)
            z_tmp = int64(CELL_DEP / 2);
        else
            z_tmp = (1 + SURFACE_LAYER);
        end
        
        tubule_cru_coordiante = tubule_cru_coordinate((tubule_cru_coordinate(:, 3) == z_tmp), :);
        for cru_id = 1 : (length(tubule_cru_coordiante(:, 1)) - 1)
            if ( (tubule_cru_coordiante(cru_id, 1)==tubule_cru_coordiante(cru_id+1, 1)) + (tubule_cru_coordiante(cru_id, 2)==tubule_cru_coordiante(cru_id+1, 2)) == 1)
                validate_tubule_at_branch(1, branch_at_no_tmp)  = validate_tubule_at_branch(1, branch_at_no_tmp) + (tubule_cru_coordiante(cru_id, 1)~=tubule_cru_coordiante(cru_id+1, 1)) * CRU_LEN;
                validate_tubule_tt_branch(1, branch_tt_no_tmp)  = validate_tubule_tt_branch(1, branch_tt_no_tmp) + (tubule_cru_coordiante(cru_id, 2)~=tubule_cru_coordiante(cru_id+1, 2)) * CRU_WID;
            else
                branch_at_no_tmp    = branch_at_no_tmp + (validate_tubule_at_branch(1, branch_at_no_tmp) ~= 0);
                branch_tt_no_tmp    = branch_tt_no_tmp + (validate_tubule_tt_branch(1, branch_tt_no_tmp) ~= 0);
            end
        end        
        branch_at_no_tmp    = branch_at_no_tmp + 1;
        branch_tt_no_tmp    = branch_tt_no_tmp + 1;

        for cru_id = 1 : (length(tubule_cru_coordiante(:, 1)) - 1)
            if(tubule_cru_coordiante(cru_id, 1)>=(1 + SURFACE_LAYER) && tubule_cru_coordiante(cru_id, 1)<=(CELL_LEN - SURFACE_LAYER))
                if(tubule_cru_coordiante(cru_id, 2)>=(1 + SURFACE_LAYER) && tubule_cru_coordiante(cru_id, 2)<=(CELL_WID - SURFACE_LAYER))
                    cru_z_branch        = cru_z_branch + 1;
                end
            end
        end
    end
    cru_x_branch        = sum(validate_tubule_at_branch) / CRU_LEN;
    cru_y_branch        = sum(validate_tubule_tt_branch) / CRU_WID;
    cru_z_branch        = int32(cru_z_branch - cru_x_branch - cru_y_branch);
    validate_cru_tubule_in_plane        = zeros(1, 3);
    
    
    if (flag_palne == 0)
        validate_cru_tubule_in_plane(1, 1)  = cru_x_branch / CELL_DEP;
        validate_cru_tubule_in_plane(1, 2)  = cru_y_branch / CELL_DEP;
        validate_cru_tubule_in_plane(1, 3)  = cru_z_branch / CELL_DEP;
    else
        validate_cru_tubule_in_plane(1, 1)  = cru_x_branch;
        validate_cru_tubule_in_plane(1, 2)  = cru_y_branch;
        validate_cru_tubule_in_plane(1, 3)  = cru_z_branch;
    end
    
    validate_tubule_z_tt_branch_tmp     = 0.3 .* ones(1, cru_z_branch);
    validate_tubule_z_tt_branch         = [validate_tubule_z_tt_branch_tmp, zeros(1, CELL_LEN * CELL_WID * CELL_DEP - length(validate_tubule_z_tt_branch_tmp))];
end

function [validate_tubule_t_index_inner, validate_tubule_t_index_whole, validate_tubule_density_inner, validate_tubule_density_whole]    = validate_tubule_calculate_density_t_index(flag_palne, validate_tubule_at_length, validate_tubule_tt_length, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN)
    inner_len       = CELL_LEN - 2 * SURFACE_LAYER;
    inner_wid       = CELL_WID - 2 * SURFACE_LAYER;
    inner_dep       = CELL_DEP - 2 * SURFACE_LAYER;
    area_xy_inner_plate     = CRU_LEN * inner_len * CRU_WID * inner_wid;
    area_inner_plates       = area_xy_inner_plate * inner_dep;
    area_xy_whole_plate     = CRU_LEN * CELL_LEN * CRU_WID * CELL_WID;
    area_whole_plates       = area_xy_whole_plate * CELL_DEP;
    area_tubule             = EXP_AT_DIAMETER_MEAN .* validate_tubule_at_length + EXP_TT_DIAMETER_MEAN .* validate_tubule_tt_length;
    
    if(flag_palne == 0)
        validate_tubule_t_index_inner   = area_tubule / area_inner_plates;
        validate_tubule_t_index_whole   = area_tubule / area_whole_plates;
        validate_tubule_density_inner   = (validate_tubule_at_length + validate_tubule_tt_length) ./ area_inner_plates;
        validate_tubule_density_whole   = (validate_tubule_at_length + validate_tubule_tt_length) ./ area_whole_plates;
    else    
        validate_tubule_t_index_inner   = area_tubule / area_xy_inner_plate;
        validate_tubule_t_index_whole   = area_tubule / area_xy_whole_plate;
        validate_tubule_density_inner   = (validate_tubule_at_length + validate_tubule_tt_length) ./ area_xy_inner_plate;
        validate_tubule_density_whole   = (validate_tubule_at_length + validate_tubule_tt_length) ./ area_xy_whole_plate;
    end
end

% To search and record all the branches in a plate
% Then to record and return the nearest distances between branches
function [validate_tubule_dist_means, validate_tubule_dist_stds, valdiate_tubule_dist_medians]   = validate_tubule_collect_inter_distance(id_file, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER)
    % To search the branch in a certain plate
    load(['TRIAL/simulation_log_' num2str(id_file) '.mat']);
    % (cru_coordinate, cru_no in single branch, branch_no, z-plate_no)
    validate_tubule_branch_coordinate   = zeros(3, CELL_LEN * CELL_WID, CELL_LEN * CELL_WID, CELL_DEP); 
%     for z_tmp = (1 + SURFACE_LAYER) : (CELL_DEP - SURFACE_LAYER)
    z_tmp = int64(CELL_DEP / 2);
%     z_tmp = (1 + SURFACE_LAYER);
        % To filter all branches (>=2CRU in a certain xy-plate)
        branch_no_tmp                   = 1;
        branch_cru_no_tmp               = 1;
        tubule_cru_coordiante_local     = tubule_cru_coordinate((tubule_cru_coordinate(:, 3) == z_tmp), :);
        for cru_id = 1 : (length(tubule_cru_coordiante_local(:, 1)) - 1)
            if ( (tubule_cru_coordiante_local(cru_id, 1)==tubule_cru_coordiante_local(cru_id+1, 1)) + (tubule_cru_coordiante_local(cru_id, 2)==tubule_cru_coordiante_local(cru_id+1, 2)) == 1)
                validate_tubule_branch_coordinate(:, branch_cru_no_tmp, branch_no_tmp, z_tmp)      = tubule_cru_coordiante_local(cru_id, :)';
                validate_tubule_branch_coordinate(:, branch_cru_no_tmp + 1, branch_no_tmp, z_tmp)  = tubule_cru_coordiante_local(cru_id+1, :)';
                branch_cru_no_tmp   = branch_cru_no_tmp + 1;
            else
                branch_cru_no_tmp   = 1;
                branch_no_tmp       = branch_no_tmp + (validate_tubule_branch_coordinate(1, 1, branch_no_tmp, z_tmp) ~= 0);
            end
        end
%     end
    
    % To search the nearest distance of each CRU in the certain branch in
    % xy-plate: validate_tubule_branch_coordinate(cru_coordinate, cru_no in single branch, branch_no, z-plate_no)
    local_meidan_dist   = zeros(1, CELL_LEN * CELL_WID * CELL_DEP);
    id_tmp              = 1;
%     for z_tmp = (1 + SURFACE_LAYER) : (CELL_DEP - SURFACE_LAYER)
    z_tmp = int64(CELL_DEP / 2);
%     z_tmp = (1 + SURFACE_LAYER);       
        for branch_no_tmp =  1 : length(find(validate_tubule_branch_coordinate(1, 1, :, z_tmp)>0))
            local_meidan_dist(1, id_tmp)    = local_search_nearest_neighbour_based_on_tubule(validate_tubule_branch_coordinate(:, :, branch_no_tmp, z_tmp), cru_state, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER);
            id_tmp  = id_tmp + 1;
        end
%     end
    validate_tubule_dist_means      = mean(local_meidan_dist(local_meidan_dist~=0));
    validate_tubule_dist_stds       = std(local_meidan_dist(local_meidan_dist~=0));
    valdiate_tubule_dist_medians    = median(local_meidan_dist(local_meidan_dist~=0));
end

% To search the nearest distance for a single branch
function local_median_dist = local_search_nearest_neighbour_based_on_tubule(validate_tubule_branch_coordinate, cru_state, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER)
    local_median_dist   = 100;
    validate_tubule_branch_coordinate_local = validate_tubule_branch_coordinate(:, (validate_tubule_branch_coordinate(1, :) > 0));
    cru_state_local     = cru_state;
    for branch_cru_no = 1: length(validate_tubule_branch_coordinate_local(1, :))
       x_tmp    = validate_tubule_branch_coordinate_local(1, branch_cru_no);
       y_tmp    = validate_tubule_branch_coordinate_local(2, branch_cru_no);
       z_tmp    = validate_tubule_branch_coordinate_local(3,  branch_cru_no);
       id_tmp   = 1 + (x_tmp - 1) + CELL_LEN * (y_tmp - 1) + CELL_LEN * CELL_WID * (z_tmp - 1);
       cru_state_local(id_tmp)  = 0;
    end
    
    inner_len       = CELL_LEN - SURFACE_LAYER;
    inner_wid       = CELL_WID - SURFACE_LAYER;
    inner_dep       = CELL_DEP - SURFACE_LAYER;
    cru_state_local = remove_leteral_tubule(SURFACE_LAYER, inner_len, inner_wid, inner_dep, cru_state_local);
    
    for branch_cru_no = 1: length(validate_tubule_branch_coordinate_local(1, :))
       x_tmp    = validate_tubule_branch_coordinate_local(1, branch_cru_no);
       y_tmp    = validate_tubule_branch_coordinate_local(2, branch_cru_no);
       z_tmp    = validate_tubule_branch_coordinate_local(3,  branch_cru_no);
       local_median_dist_tmp    = local_search_nearest_neighbour(x_tmp, y_tmp, z_tmp, cru_state_local, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID);
       local_median_dist        = min(local_median_dist, local_median_dist_tmp);
    end
end

function local_median_dist_tmp  = local_search_nearest_neighbour(x_tmp, y_tmp, z_tmp, cru_state_local, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID)
    CRU_OT                  = sqrt(CRU_LEN*CRU_LEN + CRU_WID*CRU_WID);
    local_median_dist_tmp   = 0;
    local_tt_radius         = 1;
    local_at_radius         = 1;
    local_ot_radius         = 1;
    while local_median_dist_tmp == 0
        radius_min  = min(local_tt_radius * CRU_WID, min(local_at_radius * CRU_LEN, local_ot_radius * CRU_OT) );
        if (local_tt_radius * CRU_WID) == radius_min
            % tt-step
            [cru_y1, cru_y2]        = check_y_crus_states(cru_state_local, x_tmp, y_tmp, z_tmp, local_tt_radius, CELL_LEN, CELL_WID);
            local_median_dist_tmp   = radius_min * ((cru_y1 + cru_y2) ~= 0);
            local_tt_radius         = local_tt_radius + 1;
        elseif (local_at_radius * CRU_LEN) == radius_min        
            % at-step
            [cru_x1, cru_x2]        = check_x_crus_states(cru_state_local, x_tmp, y_tmp, z_tmp, local_at_radius, CELL_LEN, CELL_WID);
            local_median_dist_tmp   = radius_min * ((cru_x1 + cru_x2) ~= 0);
            local_at_radius         = local_at_radius + 1;
        else        
            % ot-step            
            [cru_o1, cru_o2, cru_o3, cru_o4]    = check_ot_crus_states(cru_state_local, x_tmp, y_tmp, z_tmp, local_ot_radius, CELL_LEN, CELL_WID);
            local_median_dist_tmp   = radius_min * ((cru_o1 + cru_o2 + cru_o3 + cru_o4) ~= 0);
            local_ot_radius         = local_ot_radius + 1;
        end        
    end
end

% To check the status of 2 x-axis CRUs
% To check the neighbouring 2 CRUs if local_tt_radius == 1
% The outside of cell is treated as occupied
function [cru_x1, cru_x2] = check_x_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, local_at_radius, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp - local_at_radius;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    if (local_x > 0 )
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_x1      = cru_state(1, cru_id_tmp);
    else
        cru_x1  = 1;
    end
    
    local_x     = cru_x_tmp + local_at_radius;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    if (local_x < CELL_LEN )
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_x2      = cru_state(1, cru_id_tmp);
    else
        cru_x2  = 1;
    end
end

% To check the status of 2 y-axis CRUs
% To check the neighbouring 2 CRUs if local_at_radius == 1
% The outside of cell is treated as occupied
function [cru_y1, cru_y2] = check_y_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, local_tt_radius, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp - local_tt_radius;
    local_z     = cru_z_tmp;
    if (local_y > 0)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_y1      = cru_state(1, cru_id_tmp);
    else
        cru_y1  = 1;
    end
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp + local_tt_radius;
    local_z     = cru_z_tmp;
    if (local_y < CELL_WID)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_y2      = cru_state(1, cru_id_tmp);
    else
        cru_y2  = 1;
    end
end

% To check the status of 4 corner CRUs
% To check the neighbouring 4 CRUs if local_ot_radius == 1
% The outside of cell is treated as occupied
function [cru_o1, cru_o2, cru_o3, cru_o4]    = check_ot_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, local_ot_radius, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp - local_ot_radius;
    local_y     = cru_y_tmp - local_ot_radius;
    local_z     = cru_z_tmp;
    if (local_y > 0)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_o1      = cru_state(1, cru_id_tmp);
    else
        cru_o1  = 1;
    end
    local_x     = cru_x_tmp - local_ot_radius;
    local_y     = cru_y_tmp + local_ot_radius;
    local_z     = cru_z_tmp;
    if (local_y < CELL_WID)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_o2      = cru_state(1, cru_id_tmp);
    else
        cru_o2  = 1;
    end
    local_x     = cru_x_tmp + local_ot_radius;
    local_y     = cru_y_tmp - local_ot_radius;
    local_z     = cru_z_tmp;
    if (local_y < CELL_WID)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_o3      = cru_state(1, cru_id_tmp);
    else
        cru_o3  = 1;
    end
    local_x     = cru_x_tmp + local_ot_radius;
    local_y     = cru_y_tmp + local_ot_radius;
    local_z     = cru_z_tmp;
    if (local_y < CELL_WID)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_o4      = cru_state(1, cru_id_tmp);
    else
        cru_o4  = 1;
    end
end

function cru_state_local = remove_leteral_tubule(SURFACE_LAYER, inner_len, inner_wid, inner_dep, cru_state_local)
    for t_x = 1: (inner_len+2*SURFACE_LAYER)
        for t_y = 1: (inner_wid+2*SURFACE_LAYER)
            for t_z = 1: (inner_dep+2*SURFACE_LAYER)
                cru_id = 1+ (t_x-1) + (inner_len+2*SURFACE_LAYER)*(t_y-1) + (inner_len+2*SURFACE_LAYER)*(inner_wid+2*SURFACE_LAYER)*(t_z-1);
                if ( (t_x<=SURFACE_LAYER) || (t_x>(inner_len+SURFACE_LAYER)) )        cru_state_local(1, cru_id) = 0;
                elseif ( (t_y<=SURFACE_LAYER) || (t_y>(inner_wid+SURFACE_LAYER)) )    cru_state_local(1, cru_id) = 0;
                elseif ( (t_z<=SURFACE_LAYER) || (t_z>(inner_dep+SURFACE_LAYER)) )    cru_state_local(1, cru_id) = 0;
                end
            end
        end
    end
end