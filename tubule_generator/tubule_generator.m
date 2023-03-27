function tubule_generator(tubule_length_set, EXP_BRANCH_MEAN, EXP_BRANCH_STD, EXP_AT_TT, id_file, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, CRU_LEN, CRU_WID, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN)
    % Initiation of CRU dimensions and 
    % Probability of random walk direciton and beginning points    
    inner_len       = CELL_LEN - 2 * SURFACE_LAYER;
    inner_wid       = CELL_WID - 2 * SURFACE_LAYER;
    inner_dep       = CELL_DEP - 2 * SURFACE_LAYER;
    
    cru_state           = set_leteral_tubule(SURFACE_LAYER, inner_len, inner_wid, inner_dep);    

%     [probability_at, probability_tt_y, probability_tt_z]    = get_og_probability(inner_wid, inner_dep, tubule_length_set, CRU_LEN, CRU_WID, EXP_AT_TT, EXP_TT_DIAMETER_MEAN, EXP_AT_DIAMETER_MEAN);
    
    probability_at              = ( EXP_AT_TT ) / (EXP_AT_TT + 1 + 1);                          % len_AT / ( len_AT + len_TT-y + len_TT-z)
    probability_tt_y            = 1 / (EXP_AT_TT + 1 + 1);                                      % len_TT-y / ( len_AT + len_TT-y)
    probability_tt_z            = 1 / (EXP_AT_TT + 1 + 1);                                      % len_TT-z / ( len_AT + len_TT-y)
    single_cru_tubule_length    = probability_at * CRU_LEN + (1 - probability_at) * CRU_WID;    % [um / CRU]
    
    [cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp]  = find_og_tubulated_cru(probability_at, probability_tt_y, probability_tt_z, cru_state, inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID, CELL_DEP);
    cru_x_tmp       = cru_og_x_tmp;
    cru_y_tmp       = cru_og_y_tmp;
    cru_z_tmp       = cru_og_z_tmp;
    flag_first_step = 1;
    
    branch_no               = 1;
    branch_cru_num_tmp      = 0;
    branch_cru_num          = zeros(5 * inner_len * inner_wid * inner_dep, 1);
%     branch_length_um_set    = set_new_branch_length_um(single_cru_tubule_length, EXP_BRANCH_MEAN, EXP_BRANCH_STD);
%     branch_length_um_set    = EXP_BRANCH_MEAN * (EXP_AT_TT + 1 + inner_dep/inner_wid) / (EXP_AT_TT + 1) * (EXP_BRANCH_MEAN + EXP_BRANCH_MEAN / (CRU_LEN + CRU_WID) * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN) * 2 * 3) / EXP_BRANCH_MEAN;
%     branch_length_um_set    = 2 * EXP_BRANCH_MEAN * (EXP_AT_TT * (2 * EXP_TT_DIAMETER_MEAN + CRU_LEN) / CRU_LEN + 1 * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN + CRU_WID) / CRU_WID + inner_dep/inner_wid * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN + CRU_WID) / CRU_WID) / (EXP_AT_TT + 1);
    exp_converted_branch_mean   = EXP_BRANCH_MEAN * (EXP_AT_TT + 1 + inner_dep/inner_wid) / (EXP_AT_TT + 1);
    unit_branch_length_tmp      = CRU_LEN * EXP_AT_TT + CRU_WID * (1 + inner_dep/inner_wid);
    branch_length_um_set    = exp_converted_branch_mean * (exp_converted_branch_mean + exp_converted_branch_mean / unit_branch_length_tmp * 2 * (EXP_TT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) * 2 + EXP_AT_DIAMETER_MEAN * EXP_AT_TT * 2 )) / exp_converted_branch_mean;
%     branch_length_um_set    = branch_length_um_set + EXP_AT_DIAMETER_MEAN * EXP_AT_TT / (EXP_AT_TT + 1 + inner_dep/inner_wid) + EXP_TT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) / (EXP_AT_TT + 1 + inner_dep/inner_wid);
    branch_length_um_set    = branch_length_um_set + EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN;

    branch_length_um_tmp    = 0;
    branch_length_um        = zeros(inner_len * inner_wid * inner_dep, 1);
    branch_length_at_sum_um             = 0;
    branch_length_tt_y_sum_um           = 0;
    branch_length_tt_z_sum_um           = 0;
    
    tubule_cru_no           = 0;
    tubule_length_sum       = 0;    
    tubule_cru_coordinate   = zeros(5 * inner_len * inner_wid * inner_dep, 3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (tubule_length_sum < tubule_length_set)
        if check_reinitiate_branch_or_not(probability_at, probability_tt_y, probability_tt_z, branch_length_um_tmp, branch_length_um_set, cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
            % To record the former branch
            branch_cru_num(branch_no, 1)    = branch_cru_num_tmp;
            branch_length_um(branch_no, 1)  = branch_length_um_tmp;
            branch_no   = branch_no + 1;
            
            % To reset the branch length recording
            branch_cru_num_tmp      = 0;
            branch_length_um_tmp    = 0;
            branch_length_um_set    = set_new_branch_length_um(single_cru_tubule_length, EXP_AT_TT, EXP_BRANCH_MEAN, EXP_BRANCH_STD, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN, CRU_LEN, CRU_WID, inner_dep, inner_wid);

            % To start a new branch from a random proper tubulated CRU
            [cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp]  = find_og_tubulated_cru(probability_at, probability_tt_y, probability_tt_z, cru_state, inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID, CELL_DEP);
            cru_x_tmp   = cru_og_x_tmp;
            cru_y_tmp   = cru_og_y_tmp;
            cru_z_tmp   = cru_og_z_tmp;
            
            flag_first_step = 1;
            tubule_cru_no   = tubule_cru_no + 1;
            tubule_cru_coordinate(tubule_cru_no, :) = [cru_x_tmp, cru_y_tmp, cru_z_tmp]; 
            branch_cru_num_tmp                      = branch_cru_num_tmp + 1;
        end
        
        [probability_at, probability_tt_y, probability_tt_z]    = update_dynamic_probability_at(branch_length_at_sum_um, branch_length_tt_y_sum_um, branch_length_tt_z_sum_um, tubule_length_set, EXP_AT_TT, CELL_WID, CELL_DEP, CRU_LEN, CRU_WID, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN);
        [cru_id, cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step_in_myocyte(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, probability_at, probability_tt_y, probability_tt_z, CELL_LEN, CELL_WID, CELL_DEP);
        if (cru_state(1,cru_id) == 0)
            tubule_cru_no   = tubule_cru_no + 1;
            tubule_cru_coordinate(tubule_cru_no, :) = [cru_x_tmp, cru_y_tmp, cru_z_tmp];
            
            if(flag_first_step == 1)
                flag_first_step = 0;
                if ( cru_x_tmp ~= cru_og_x_tmp )
%                     branch_length_at_sum_um     = branch_length_at_sum_um + CRU_LEN;
%                     branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + 2 * EXP_TT_DIAMETER_MEAN;
%                     branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + 2 * EXP_TT_DIAMETER_MEAN;
%                     branch_length_um_tmp    = branch_length_um_tmp + CRU_LEN + 4 * EXP_TT_DIAMETER_MEAN;
%                     tubule_length_sum       = tubule_length_sum + CRU_LEN + 4 * EXP_TT_DIAMETER_MEAN;
                    branch_length_at_sum_um     = branch_length_at_sum_um + CRU_LEN;
                    branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + EXP_TT_DIAMETER_MEAN;
                    branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + EXP_TT_DIAMETER_MEAN;
                    branch_length_um_tmp    = branch_length_um_tmp + CRU_LEN + 2 * EXP_TT_DIAMETER_MEAN;
                    tubule_length_sum       = tubule_length_sum + CRU_LEN + 2 * EXP_TT_DIAMETER_MEAN;
                else
                    if ( cru_y_tmp ~= cru_og_y_tmp )
%                         branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + CRU_WID;
%                         branch_length_at_sum_um     = branch_length_at_sum_um + 2 * EXP_AT_DIAMETER_MEAN;
%                         branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + 2 * EXP_TT_DIAMETER_MEAN;
                        branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + CRU_WID;
                        branch_length_at_sum_um     = branch_length_at_sum_um + EXP_AT_DIAMETER_MEAN;
                        branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + EXP_TT_DIAMETER_MEAN;
                    else
%                         branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + CRU_WID;
%                         branch_length_at_sum_um     = branch_length_at_sum_um + 2 * EXP_AT_DIAMETER_MEAN;
%                         branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + 2 * EXP_TT_DIAMETER_MEAN;
                        branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + CRU_WID;
                        branch_length_at_sum_um     = branch_length_at_sum_um + EXP_AT_DIAMETER_MEAN;
                        branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + EXP_TT_DIAMETER_MEAN;
                    end
%                     branch_length_um_tmp        = branch_length_um_tmp + CRU_WID + 2 * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN);
%                     tubule_length_sum           = tubule_length_sum + CRU_WID + 2 * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN);
                    branch_length_um_tmp        = branch_length_um_tmp + CRU_WID + EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN;
                    tubule_length_sum           = tubule_length_sum + CRU_WID + EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN;
                end
            else    
                if ( cru_x_tmp ~= tubule_cru_coordinate(tubule_cru_no-1, 1) )
                    branch_length_at_sum_um     = branch_length_at_sum_um + CRU_LEN;
                    branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + EXP_TT_DIAMETER_MEAN;
                    branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + EXP_TT_DIAMETER_MEAN;
                    branch_length_um_tmp    = branch_length_um_tmp + CRU_LEN + 2 * EXP_TT_DIAMETER_MEAN;
                    tubule_length_sum       = tubule_length_sum + CRU_LEN + 2 * EXP_TT_DIAMETER_MEAN;
                else
                    if ( cru_y_tmp ~= tubule_cru_coordinate(tubule_cru_no-1, 2) )
                        branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + CRU_WID;
                        branch_length_at_sum_um     = branch_length_at_sum_um + EXP_AT_DIAMETER_MEAN;
                        branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + EXP_TT_DIAMETER_MEAN;
                    else
                        branch_length_tt_z_sum_um   = branch_length_tt_z_sum_um + CRU_WID;
                        branch_length_at_sum_um     = branch_length_at_sum_um + EXP_AT_DIAMETER_MEAN;
                        branch_length_tt_y_sum_um   = branch_length_tt_y_sum_um + EXP_TT_DIAMETER_MEAN;
                    end
                    branch_length_um_tmp        = branch_length_um_tmp + CRU_WID + EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN;
                    tubule_length_sum           = tubule_length_sum + CRU_WID + EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN;
                end
            end

            cru_id  = 1 + (cru_x_tmp - 1) + CELL_LEN * (cru_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_z_tmp - 1);
            cru_state(1, cru_id)    = 1;
            branch_cru_num_tmp      = branch_cru_num_tmp + 1;
        end
        
        if (tubule_cru_no > 0) && (flag_first_step ~= 1)
            cru_x_tmp   = tubule_cru_coordinate(tubule_cru_no, 1);
            cru_y_tmp   = tubule_cru_coordinate(tubule_cru_no, 2);
            cru_z_tmp   = tubule_cru_coordinate(tubule_cru_no, 3);
        elseif (flag_first_step == 1)
            cru_x_tmp   = cru_og_x_tmp;
            cru_y_tmp   = cru_og_y_tmp;
            cru_z_tmp   = cru_og_z_tmp;
        end
        
        if (sum(cru_state) == length(cru_state) )
            break
        end
    end

    branch_cru_num(branch_no, 1)    = branch_cru_num_tmp;
    branch_length_um(branch_no, 1)  = branch_length_um_tmp;
    save(['TRIAL/simulation_log_' num2str(id_file) '.mat'])
    
%     tubule_structure_visul(id_file, tubule_cru_no, tubule_cru_coordinate);
    disp(id_file);
    cru_input(cru_state, id_file);
    cru_paraview_visul(cru_state, id_file, CELL_LEN, CELL_WID, CELL_DEP, CRU_LEN, CRU_WID);
    tubule_paraview_visul(tubule_cru_coordinate, tubule_cru_no, branch_cru_num, branch_no, id_file);
end

function flag = check_reinitiate_branch_or_not(probability_at, probability_tt_y, probability_tt_z, branch_length_um_tmp, branch_length_um_set, cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
    flag1   = (branch_length_um_tmp >= branch_length_um_set);
    flag2   = check_around_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP);
    flag3   = ( probability_at == 1 && check_around_at_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID) );
    flag4   = ( probability_at == 0 && check_around_tt_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) );
    flag5   = ( probability_tt_y == 1 && check_around_tt_y_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID) );
    flag6_1 = ( probability_tt_y == 0 );
    flag6_2 = ( check_around_at_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID) && check_around_tt_z_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) );
    flag6   = flag6_1 && flag6_2;
    flag7   = ( probability_tt_z == 1 && check_around_tt_z_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) );
    flag8_1 = ( probability_tt_z == 0 );
    flag8_2 = ( check_around_at_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID) && check_around_tt_y_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID) );
    flag8   = ( flag8_1 && flag8_2 );
    
    flag    = ( flag1 || flag2 || flag3 || flag4 || flag5 || flag6 || flag7 || flag8 );
end

% To set up a new branch length and rescale from 2D branch length to 3D
% by adding potential z-direnction branch length
function branch_length_um_set = set_new_branch_length_um(single_cru_tubule_length, EXP_AT_TT, EXP_BRANCH_MEAN, EXP_BRANCH_STD, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN, CRU_LEN, CRU_WID, inner_dep, inner_wid)
%     branch_length_um_set = 0;
%     while ( (branch_length_um_set <= single_cru_tubule_length) || (branch_length_um_set >= (2 * EXP_BRANCH_MEAN - single_cru_tubule_length) ) )
        branch_length_um_set = normrnd(EXP_BRANCH_MEAN, EXP_BRANCH_STD);
%     end

    if (branch_length_um_set <= single_cru_tubule_length)
        branch_length_um_set = single_cru_tubule_length;
%     elseif (branch_length_um_set >= (2 * EXP_BRANCH_MEAN - single_cru_tubule_length) )
%         branch_length_um_set = 2 * EXP_BRANCH_MEAN - single_cru_tubule_length;
    end
    
    branch_length_um_set    = branch_length_um_set * (EXP_AT_TT + 1 + inner_dep/inner_wid) / (EXP_AT_TT + 1);
    exp_converted_branch_mean   = EXP_BRANCH_MEAN * (EXP_AT_TT + 1 + inner_dep/inner_wid) / (EXP_AT_TT + 1);
    unit_branch_length_tmp      = CRU_LEN * EXP_AT_TT + CRU_WID * (1 + inner_dep/inner_wid);
%     branch_length_um_set    = branch_length_um_set * (1.0 + exp_converted_branch_mean / unit_branch_length_tmp * 2 * ( EXP_AT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) * 2 + EXP_TT_DIAMETER_MEAN * EXP_AT_TT * 2 ));
    
%     branch_length_um_set    = branch_length_um_set * (1.0 + 1.0 / unit_branch_length_tmp * 2 * ( EXP_TT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) * 2 + EXP_AT_DIAMETER_MEAN * EXP_AT_TT * 2 ));
    branch_length_um_set    = branch_length_um_set * (1.0 + 1.0 / unit_branch_length_tmp * 2 * ( EXP_TT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) * 2 + EXP_AT_DIAMETER_MEAN * EXP_AT_TT * 2 ));
%     branch_length_um_set    = branch_length_um_set + EXP_AT_DIAMETER_MEAN * EXP_AT_TT / (EXP_AT_TT + 1 + inner_dep/inner_wid) + EXP_TT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) / (EXP_AT_TT + 1 + inner_dep/inner_wid);
    branch_length_um_set    = branch_length_um_set + EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN;
    
%     branch_length_um_set =2 * branch_length_um_set * (EXP_AT_TT * (2 * EXP_TT_DIAMETER_MEAN + CRU_LEN) / CRU_LEN + 1 * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN + CRU_WID) / CRU_WID + inner_dep/inner_wid * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN + CRU_WID) / CRU_WID) / (EXP_AT_TT + 1);
    
%     if (branch_length_um_set <= 0)
%         branch_length_um_set = 0;
%     elseif (branch_length_um_set >= (2 * EXP_BRANCH_MEAN) )
%         branch_length_um_set = 2 * EXP_BRANCH_MEAN;
%     end
end

% To set the surface 2-layer CRUs coupled with tubules
% CRU is ordered in: cru_id = cru_x_tmp + inner_len * cru_y_tmp + inner_len * inner_wid * cru_z_tmp;
% State: Coupled - 1; uncoupled - 0
function cru_state = set_leteral_tubule(SURFACE_LAYER, inner_len, inner_wid, inner_dep)
    cru_state     = zeros(1, (inner_len+2*SURFACE_LAYER)*(inner_wid+2*SURFACE_LAYER)*(inner_dep+2*SURFACE_LAYER) );    
    for t_x = 1: (inner_len+2*SURFACE_LAYER)
        for t_y = 1: (inner_wid+2*SURFACE_LAYER)
            for t_z = 1: (inner_dep+2*SURFACE_LAYER)
                cru_id = 1+ (t_x-1) + (inner_len+2*SURFACE_LAYER)*(t_y-1) + (inner_len+2*SURFACE_LAYER)*(inner_wid+2*SURFACE_LAYER)*(t_z-1);
                if ( (t_x<=SURFACE_LAYER) || (t_x>(inner_len+SURFACE_LAYER)) )        cru_state(1, cru_id) = 1;
                elseif ( (t_y<=SURFACE_LAYER) || (t_y>(inner_wid+SURFACE_LAYER)) )    cru_state(1, cru_id) = 1;
                elseif ( (t_z<=SURFACE_LAYER) || (t_z>(inner_dep+SURFACE_LAYER)) )    cru_state(1, cru_id) = 1;
                end
            end
        end
    end
end

% To find a tubulated CRU with un-tubulated enighbouring CRU
% This tubuldated CRU can be a inner-layer CRU
function [cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = find_og_tubulated_cru(probability_at, probability_tt_y, probability_tt_z, cru_state, inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID, CELL_DEP)
    [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
    if (probability_at > 0) && (probability_at < 1) && (probability_tt_y > 0) && (probability_tt_y < 1) && (probability_tt_z > 0) && (probability_tt_z < 1)
        while ( (cru_state(1, cru_og_id) == 0) || check_around_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    elseif (probability_at == 1)
        while ( (cru_state(1, cru_og_id) == 0) || check_around_at_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    elseif (probability_tt_y == 1)
        while ( (cru_state(1, cru_og_id) == 0) || check_around_tt_y_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    elseif (probability_tt_z == 1)
        while ( (cru_state(1, cru_og_id) == 0) || check_around_tt_z_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    elseif (probability_at == 0)
        while ( (cru_state(1, cru_og_id) == 0) || check_around_tt_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    elseif (probability_tt_y == 0)
        while ( (cru_state(1, cru_og_id) == 0) || (check_around_at_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID) && check_around_tt_z_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) ) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    elseif (probability_tt_z == 0)
        while ( (cru_state(1, cru_og_id) == 0) || (check_around_at_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID) && check_around_tt_y_crus_all_occupied(cru_state, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, CELL_LEN, CELL_WID) ) )
            [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID);
        end
    end
end

function [cru_og_id, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp] = random_select_cru(inner_len, inner_wid, inner_dep, SURFACE_LAYER, CELL_LEN, CELL_WID)
    cru_og_x_tmp    = randi(inner_len + 2) + SURFACE_LAYER - 1;
    cru_og_y_tmp    = randi(inner_wid + 2) + SURFACE_LAYER - 1;
    cru_og_z_tmp    = randi(inner_dep + 2) + SURFACE_LAYER - 1;
    cru_og_id       = 1 + (cru_og_x_tmp - 1) + CELL_LEN * (cru_og_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_og_z_tmp - 1);
end

% To check the status of centeral CRU and all 2 AT CRUs around
function flag = check_around_at_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
    cru_og  = cru_state(1, cru_id_tmp);
    
    [cru_x1, cru_x2] = check_x_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID);
    flag    = ( (cru_og * cru_x1 * cru_x2) == 1 );
end

% To check the status of centeral CRU and all 2 y-TT CRUs around
function flag = check_around_tt_y_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
    cru_og  = cru_state(1, cru_id_tmp);
    
    [cru_y1, cru_y2] = check_y_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID);
    flag    = ( (cru_og * cru_y1 * cru_y2) == 1 );
end

% To check the status of centeral CRU and all 2 z-TT CRUs around
function flag = check_around_tt_z_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
    cru_og  = cru_state(1, cru_id_tmp);
    
    [cru_z1, cru_z2] = check_z_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP);
    flag    = ( (cru_og * cru_z1 * cru_z2) == 1 );
end

% To check the status of centeral CRU and all 4 TT CRUs around
% If all the 7 CRUs are occupied, new branch gets initiated
function flag = check_around_tt_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
    flag    = ( check_around_tt_y_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID) * check_around_tt_z_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP) );
end

% To check the status of centeral CRU and all 6 CRUs around
% If all the 7 CRUs are occupied, new branch gets initiated
function flag = check_around_crus_all_occupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
    cru_og  = cru_state(1, cru_id_tmp);
    
    [cru_x1, cru_x2, cru_y1, cru_y2, cru_z1, cru_z2]    = check_around_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP);
    flag    = ( (cru_og * cru_x1 * cru_x2 * cru_y1 * cru_y2 * cru_z1 * cru_z2) == 1 );
end

% % To check the status of all 6 CRUs around
% % If all the 6 CRUs are unoccupied, new branch is re-initiated
% function flag = check_around_crus_all_unoccupied(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
%     [cru_x1, cru_x2, cru_y1, cru_y2, cru_z1, cru_z2]    = check_around_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP);
%     flag    = ( (cru_x1 + cru_x2 + cru_y1 + cru_y2 + cru_z1 + cru_z2) == 0 );
% end

% To check the status of 2 neighbouring x-axis CRUs around
% The outside of cell is treated as occupied
function [cru_x1, cru_x2] = check_x_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp - 1;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    if (local_x > 0 )
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_x1      = cru_state(1, cru_id_tmp);
    else
        cru_x1  = 1;
    end
    
    local_x     = cru_x_tmp + 1;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp;
    if (local_x < CELL_LEN )
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_x2      = cru_state(1, cru_id_tmp);
    else
        cru_x2  = 1;
    end
end

% To check the status of all 2 neighbouring y-CRUs around
% The outside of cell is treated as occupied
function [cru_y1, cru_y2] = check_y_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp - 1;
    local_z     = cru_z_tmp;
    if (local_y > 0)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_y1      = cru_state(1, cru_id_tmp);
    else
        cru_y1  = 1;
    end
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp + 1;
    local_z     = cru_z_tmp;
    if (local_y < CELL_WID)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_y2      = cru_state(1, cru_id_tmp);
    else
        cru_y2  = 1;
    end
end

% To check the status of all 2 neighbouring z-CRUs around
% The outside of cell is treated as occupied
function [cru_z1, cru_z2] = check_z_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp - 1;
    if (local_z > 0)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_z1      = cru_state(1, cru_id_tmp);
    else
        cru_z1  = 1;
    end
    
    local_x     = cru_x_tmp;
    local_y     = cru_y_tmp;
    local_z     = cru_z_tmp + 1;
    if (local_z < CELL_DEP)
        cru_id_tmp  = 1 + (local_x - 1) + CELL_LEN * (local_y - 1) + CELL_LEN * CELL_WID * (local_z - 1);
        cru_z2      = cru_state(1, cru_id_tmp);
    else
        cru_z2  = 1;
    end
end

% To check the status of all 6 CRUs around
% The outside of cell is treated as occupied
function [cru_x1, cru_x2, cru_y1, cru_y2, cru_z1, cru_z2] = check_around_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP)
    [cru_x1, cru_x2] = check_x_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID);
    [cru_y1, cru_y2] = check_y_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID);
    [cru_z1, cru_z2] = check_z_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP);
end

function [probability_at, probability_tt_y, probability_tt_z] = update_dynamic_probability_at(branch_length_at_sum_um, branch_length_tt_y_sum_um, branch_length_tt_z_sum_um, tubule_length_set, EXP_AT_TT, CELL_WID, CELL_DEP, CRU_LEN, CRU_WID, EXP_AT_DIAMETER, EXP_TT_DIAMETER)
    probability_at_og   = ( EXP_AT_TT ) / (EXP_AT_TT + 1 + 1);          % #AT / ( #AT + #TT-y)
    probability_tt_y_og = 1 / (EXP_AT_TT + 1 + 1);                      % len_TT-y / ( len_AT + len_TT-y)
    probability_tt_z_og = 1 / (EXP_AT_TT + 1 + 1);                      % len_TT-z / ( len_AT + len_TT-y)
%     [probability_at_og, probability_tt_y_og, probability_tt_z_og]    = get_og_probability(inner_wid, inner_dep, tubule_length_set, CRU_LEN, CRU_WID, EXP_AT_TT, EXP_TT_DIAMETER_MEAN, EXP_AT_DIAMETER_MEAN);
    remain_at   = (probability_at_og * tubule_length_set - branch_length_at_sum_um);
    remain_tt_y = (probability_tt_y_og * tubule_length_set - branch_length_tt_y_sum_um);
    remain_tt_z = (probability_tt_z_og * tubule_length_set - branch_length_tt_z_sum_um);
    if ( (branch_length_tt_y_sum_um + branch_length_tt_z_sum_um) >= tubule_length_set * (1 - probability_at_og) )
        probability_at      = 1;
        probability_tt_y    = 0;
        probability_tt_z    = 0;
    elseif ( (branch_length_at_sum_um + branch_length_tt_z_sum_um) >= tubule_length_set * (1 - probability_tt_y_og) )
        probability_at      = 0;
        probability_tt_y    = 1;
        probability_tt_z    = 0;
    elseif ( (branch_length_at_sum_um + branch_length_tt_y_sum_um) >= tubule_length_set * (1 - probability_tt_z_og) )
        probability_at      = 0;
        probability_tt_y    = 0;
        probability_tt_z    = 1;
    elseif (branch_length_at_sum_um >= tubule_length_set * probability_at_og)
        probability_at      = 0;
        probability_tt_y    = (probability_tt_y_og * tubule_length_set - branch_length_tt_y_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
        probability_tt_z    = (probability_tt_z_og * tubule_length_set - branch_length_tt_z_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
    elseif (branch_length_tt_y_sum_um >= tubule_length_set * probability_tt_y_og)
        probability_at      = (probability_at_og * tubule_length_set - branch_length_at_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
        probability_tt_y    = 0;
        probability_tt_z    = (probability_tt_z_og * tubule_length_set - branch_length_tt_z_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
    elseif (branch_length_tt_z_sum_um >= tubule_length_set * probability_tt_z_og)
        probability_at      = (probability_at_og * tubule_length_set - branch_length_at_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
        probability_tt_y    = (probability_tt_y_og * tubule_length_set - branch_length_tt_y_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
        probability_tt_z    = 0;
    else
        probability_at      = (probability_at_og * tubule_length_set - branch_length_at_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
        probability_tt_y    = (probability_tt_y_og * tubule_length_set - branch_length_tt_y_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
        probability_tt_z    = (probability_tt_z_og * tubule_length_set - branch_length_tt_z_sum_um) / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um );
    end
%      && remain_tt_y / CRU_WID * EXP_TT_DIAMETER < remain_tt_z && remain_tt_y / CRU_WID * EXP_TT_DIAMETER < remain_at 
%  && remain_tt_z / CRU_WID * EXP_TT_DIAMETER < remain_tt_y && remain_tt_z / CRU_WID * EXP_TT_DIAMETER < remain_at
%     if ( remain_at / CRU_LEN * EXP_AT_DIAMETER < remain_tt_y || remain_at / CRU_LEN * EXP_AT_DIAMETER < remain_tt_z )
%         probability_at      = 0;
%         probability_tt_y    = abs(remain_tt_y / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um ));
%         probability_tt_z    = abs(remain_tt_z / ( tubule_length_set - branch_length_at_sum_um - branch_length_tt_y_sum_um - branch_length_tt_z_sum_um ));
%     end
        
end

% To grow one step and make sure tubule will not get out of the whole cell
function [cru_id, cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step_in_myocyte(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, cru_og_x_tmp, cru_og_y_tmp, cru_og_z_tmp, probability_at, probability_tt_y, probability_tt_z, CELL_LEN, CELL_WID, CELL_DEP)
    [cru_x1, cru_x2, cru_y1, cru_y2, cru_z1, cru_z2] = check_around_crus_states(cru_state, cru_x_tmp, cru_y_tmp, cru_z_tmp, CELL_LEN, CELL_WID, CELL_DEP);
    if (cru_x1 + cru_x2 + cru_y1 + cru_y2 + cru_z1 + cru_z2) == 5
        cru_x_tmp   = cru_x_tmp - (cru_x1 == 0) + (cru_x2 == 0);
        cru_y_tmp   = cru_y_tmp - (cru_y1 == 0) + (cru_y2 == 0);
        cru_z_tmp   = cru_z_tmp - (cru_z1 == 0) + (cru_z2 == 0);
        cru_id      = 1 + (cru_x_tmp - 1) + CELL_LEN * (cru_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_z_tmp - 1);
    elseif (cru_x1 + cru_x2) == 2
        [cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step_yzplate(cru_x_tmp, cru_y_tmp, cru_z_tmp, probability_tt_y, probability_tt_z);
        cru_id  = 1 + (cru_x_tmp - 1) + CELL_LEN * (cru_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_z_tmp - 1);
        while (cru_x_tmp < 1) || (cru_x_tmp > CELL_LEN) || (cru_y_tmp < 1) || (cru_y_tmp > CELL_WID) || (cru_z_tmp < 1) || (cru_z_tmp > CELL_DEP)
            cru_x_tmp   = cru_og_x_tmp;
            cru_y_tmp   = cru_og_y_tmp;
            cru_z_tmp   = cru_og_z_tmp;
            [cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step(cru_x_tmp, cru_y_tmp, cru_z_tmp, probability_at, probability_tt_y, probability_tt_z);
            cru_id  = 1 + (cru_x_tmp - 1) + CELL_LEN * (cru_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_z_tmp - 1);
        end
    else
        [cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step(cru_x_tmp, cru_y_tmp, cru_z_tmp, probability_at, probability_tt_y, probability_tt_z);
        cru_id  = 1 + (cru_x_tmp - 1) + CELL_LEN * (cru_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_z_tmp - 1);
        while (cru_x_tmp < 1) || (cru_x_tmp > CELL_LEN) || (cru_y_tmp < 1) || (cru_y_tmp > CELL_WID) || (cru_z_tmp < 1) || (cru_z_tmp > CELL_DEP)
            cru_x_tmp   = cru_og_x_tmp;
            cru_y_tmp   = cru_og_y_tmp;
            cru_z_tmp   = cru_og_z_tmp;
            [cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step(cru_x_tmp, cru_y_tmp, cru_z_tmp, probability_at, probability_tt_y, probability_tt_z);
            cru_id  = 1 + (cru_x_tmp - 1) + CELL_LEN * (cru_y_tmp - 1) + CELL_LEN * CELL_WID * (cru_z_tmp - 1);
        end
    end
end

function [cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step_yzplate(cru_x_tmp, cru_y_tmp, cru_z_tmp, probability_tt_y, probability_tt_z)
        r   = rand(1);
        if ( r < probability_tt_y/2 / (probability_tt_y+probability_tt_z) )
            cru_y_tmp  = cru_y_tmp - 1;
        elseif ( r < probability_tt_y / (probability_tt_y+probability_tt_z) )
            cru_y_tmp  = cru_y_tmp + 1;
        elseif ( r < (probability_tt_y+probability_tt_z/2) / (probability_tt_y+probability_tt_z) )
            cru_z_tmp  = cru_z_tmp - 1;
        else
            cru_z_tmp  = cru_z_tmp + 1;
        end
end

function [cru_x_tmp, cru_y_tmp, cru_z_tmp]   = tubule_random_grow_one_step(cru_x_tmp, cru_y_tmp, cru_z_tmp, probability_at, probability_tt_y, probability_tt_z)
        r   = rand(1);
        if ( r < probability_at/2 )
            cru_x_tmp  = cru_x_tmp - 1;
        elseif ( r < probability_at )
            cru_x_tmp  = cru_x_tmp + 1;
        elseif ( r < probability_at + probability_tt_y/2 )
            cru_y_tmp  = cru_y_tmp - 1;
        elseif ( r < probability_at + probability_tt_y )
            cru_y_tmp  = cru_y_tmp + 1;
        elseif ( r < probability_at + probability_tt_y + probability_tt_z/2 )
            cru_z_tmp  = cru_z_tmp - 1;
        else
            cru_z_tmp  = cru_z_tmp + 1;
        end
end

function tubule_structure_visul(id_file, tubule_cru_no, tubule_cru_coordinate)
    figure(id_file)
    scatter3(tubule_cru_coordinate(1:tubule_cru_no,1), tubule_cru_coordinate(1:tubule_cru_no, 2), tubule_cru_coordinate(1:tubule_cru_no, 3), 150, 'filled')
    axis equal
end

function cru_paraview_visul(cru_state, id_file, CELL_LEN, CELL_WID, CELL_DEP, CRU_LEN, CRU_WID)
    n_CRU   = CELL_LEN*CELL_WID*CELL_DEP;
    filename_cru = fopen(['TRIAL/cru_visul_' num2str(id_file) '.vtk'],'w');
    fprintf(filename_cru, '# vtk DataFile Version 3.0\n');
    fprintf(filename_cru, 'vtk output\n');
    fprintf(filename_cru, 'ASCII\n');
    fprintf(filename_cru, 'DATASET STRUCTURED_POINTS\n');
    fprintf(filename_cru, 'DIMENSIONS %d %d %d\n', CELL_LEN, CELL_WID, CELL_DEP);
    fprintf(filename_cru, 'SPACING %.2f %.1f %.1f\n', CRU_LEN, CRU_WID, CRU_WID);
    fprintf(filename_cru, 'ORIGIN 0 0 0\n');
    fprintf(filename_cru, 'POINT_DATA %d\n', n_CRU);
    fprintf(filename_cru, 'SCALARS HumanAtrium float 1\n');
    fprintf(filename_cru, 'LOOKUP_TABLE default\n');
    fprintf(filename_cru, '%5d', cru_state);
    fclose(filename_cru);
end

% To output the cru state as simulation input
function cru_input(cru_state, id_file)
    filename_cru_input  = fopen(['TRIAL/tub_input_ver2_' num2str(id_file) '.txt'], 'w');
    for local_cru_no = 1:length(cru_state)
        if local_cru_no < length(cru_state)
            fprintf(filename_cru_input, '%d\n\r', cru_state(1, local_cru_no) );
        else
            fprintf(filename_cru_input, '%d', cru_state(1, local_cru_no) );
        end 
    end
    fclose(filename_cru_input);
end

% To output the tubule sturcture using specific format in paraview
function tubule_paraview_visul(tubule_cru_coordinate, tubule_cru_no, branch_cru_num, branch_no, id_file)    
    filename_tubule_input   = fopen(['TRIAL/tubule_visul_' num2str(id_file) '.vtk'], 'w');
    fprintf(filename_tubule_input, '# vtk DataFile Version 2.0\n');
    fprintf(filename_tubule_input, 'Unstructured Grid Example\n');
    fprintf(filename_tubule_input, 'ASCII\n');
    fprintf(filename_tubule_input, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(filename_tubule_input, '\n');
    fprintf(filename_tubule_input, 'POINTS %d float\n', tubule_cru_no);

%     for i = 1 : length(cru_state)
%        if ( (branch_cru_num(i, 1) > 0) && (branch_cru_num(i + 1, 1) < 1) )
%            branch_num   = i;
%            break
%        end
%     end
    

    tubule_cru_coordinate_local     = tubule_cru_coordinate - 1;
    tubule_cru_no_local             = 0;
%     for i = 1 : branch_num
    for i = 1 : branch_no
       for j = 1 : branch_cru_num(i, 1)
            tubule_cru_no_local     = tubule_cru_no_local + 1;
            fprintf(filename_tubule_input, '%d %d %d  ', tubule_cru_coordinate_local(tubule_cru_no_local, 1), tubule_cru_coordinate_local(tubule_cru_no_local, 2), tubule_cru_coordinate_local(tubule_cru_no_local, 3) );   
       end
       fprintf(filename_tubule_input, '\n');
    end

    fprintf(filename_tubule_input, '\n');
    fprintf(filename_tubule_input, 'CELLS %d %d\n', branch_no, (branch_no + tubule_cru_no ));
    cell_order_local    = 0;
    for i = 1 : branch_no
        fprintf( filename_tubule_input, '%d ', branch_cru_num(i, 1) );
        for j = 1:branch_cru_num(i)
            fprintf(filename_tubule_input, '%d ', cell_order_local );
            cell_order_local    = cell_order_local + 1;
        end
        fprintf(filename_tubule_input, '\n');
    end

    fprintf(filename_tubule_input, '\n');
    fprintf(filename_tubule_input, 'CELL_TYPES %d\n', branch_no );
    for i = 1 : branch_no
        if (branch_cru_num(i, 1) == 1)
            fprintf(filename_tubule_input, '%d\n', 1);
        else
            fprintf(filename_tubule_input, '%d\n', 3);
        end
    end

    fprintf(filename_tubule_input, '\n');
    fprintf(filename_tubule_input, 'POINT_DATA %d\n', tubule_cru_no );
    fprintf(filename_tubule_input, 'SCALARS scalars float 1\n');
    fprintf(filename_tubule_input, 'LOOKUP_TABLE default\n');
    for i = 1 : tubule_cru_no
        fprintf(filename_tubule_input, '1.0 ');
        if ( mode(i, 10) == 0 )
            fprintf(filename_tubule_input, '\n'); 
        end
    end
    fclose(filename_tubule_input);
end

function [probability_at_og, probability_tt_y_og, probability_tt_z_og]    = get_og_probability(inner_wid, inner_dep, tubule_length_set, CRU_LEN, CRU_WID, EXP_AT_TT, EXP_TT_DIAMETER_MEAN, EXP_AT_DIAMETER_MEAN)
    tubule_length_x     = ( EXP_AT_TT ) / (EXP_AT_TT + 1 + 1) * tubule_length_set;
    tubule_length_y     = 1.0 / (EXP_AT_TT + 1 + 1) * tubule_length_set;
    tubule_length_z     = tubule_length_y;
    
    syms AT_num TT_y_num TT_z_num
%     equ1    = tubule_length_x - AT_num * CRU_LEN - (TT_y_num + TT_z_num) * (EXP_TT_DIAMETER_MEAN * 2);
%     equ2    = tubule_length_y - AT_num * (EXP_AT_DIAMETER_MEAN * 2) - TT_y_num * CRU_WID - TT_z_num * (EXP_TT_DIAMETER_MEAN * 2);
%     equ3    = tubule_length_z - AT_num * (EXP_AT_DIAMETER_MEAN * 2) - TT_y_num * (EXP_TT_DIAMETER_MEAN * 2) - TT_z_num * CRU_WID;
%     [AT_num, TT_y_num, TT_z_num]      = solve(equ1, equ2, equ3);
    
    Y = vpasolve([0 == tubule_length_x - AT_num * CRU_LEN - (TT_y_num + TT_z_num) * (EXP_TT_DIAMETER_MEAN), 0 == tubule_length_y - AT_num * (EXP_AT_DIAMETER_MEAN) - TT_y_num * CRU_WID - TT_z_num * (EXP_TT_DIAMETER_MEAN), 0 == tubule_length_z - AT_num * (EXP_AT_DIAMETER_MEAN) - TT_y_num * (EXP_TT_DIAMETER_MEAN) - TT_z_num * CRU_WID], [AT_num, TT_y_num, TT_z_num]);
    
    probability_at_og   = Y.AT_num / (Y.AT_num + Y.TT_y_num + Y.TT_z_num);
    probability_tt_y_og = TT_y_num / (AT_num + TT_y_num + TT_z_num);
    probability_tt_z_og = TT_z_num / (AT_num + TT_y_num + TT_z_num);
end