% Tubular pattern generator
% Use random walk to generate a pupultation of TAT strucutres
% Using experimental data (tubular density, branch length, ratio pf AT/TT)
% Update 04/10/2020, Xianwei Zhang <xianweizhang95@gmail.com>
function tubule_population_generator()   
    clear
    clc
    
    EXP_DEN_MEAN    = 0.37;
    EXP_DEN_STD     = 0.17;
    EXP_AT_TT       = 6;
    EXP_BRANCH_MEAN = 2.29;
    EXP_BRANCH_STD  = 0.98;
    EXP_AT_DIAMETER_MEAN    = 0.32426;
    EXP_TT_DIAMETER_MEAN    = 0.27219;

    CELL_LEN        = 55;
    CELL_WID        = 17;
    CELL_DEP        = 11;
    SURFACE_LAYER   = 2;
    CRU_LEN         = 1.84;
    CRU_WID         = 0.9;
    
    TUBULE_POPULATION_NUM   = 19;
    tubule_length_set   = set_tubule_length_um(TUBULE_POPULATION_NUM, SURFACE_LAYER, EXP_DEN_MEAN, EXP_DEN_STD, EXP_BRANCH_MEAN, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN, CRU_LEN, CRU_WID, EXP_AT_TT, CELL_LEN, CELL_WID, CELL_DEP);
    for id_file = 1 : (TUBULE_POPULATION_NUM)
        tubule_generator(tubule_length_set(id_file, 1), EXP_BRANCH_MEAN, EXP_BRANCH_STD, EXP_AT_TT, id_file, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, CRU_LEN, CRU_WID, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN);
    end
    disp('Building population Done!');

    validate_based_on_tubule_state_single_plane_involve_z_branch(TUBULE_POPULATION_NUM, CRU_LEN, CRU_WID, CELL_LEN, CELL_WID, CELL_DEP, SURFACE_LAYER, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN);
    disp('Tubule structure - single plane - involve z-branches - Validation Done!');
    
    show_exp_sim_comparison(CELL_LEN, CELL_WID, CELL_DEP, TUBULE_POPULATION_NUM);
    disp('Figure Plotting - Done!');
end

% Use experimental data (tubular density, size of cell and CRU) to
% calculate and set the whole length [um] of tubule in single cell
% and add z-tubule (2D -> 3D density)
function tubule_length_set = set_tubule_length_um(TUBULE_POPULATION_NUM, SURFACE_LAYER, EXP_DEN_MEAN, EXP_DEN_STD, EXP_BRANCH_MEAN, EXP_AT_DIAMETER_MEAN, EXP_TT_DIAMETER_MEAN, CRU_LEN, CRU_WID, EXP_AT_TT, CELL_LEN, CELL_WID, CELL_DEP)
    tubule_length_set   = zeros(TUBULE_POPULATION_NUM + 1, 1);
    inner_len       = CELL_LEN - 2 * SURFACE_LAYER;
    inner_wid       = CELL_WID - 2 * SURFACE_LAYER;
    inner_dep       = CELL_DEP - 2 * SURFACE_LAYER;
    tubule_density_per_cru_2d_mean  = EXP_DEN_MEAN * CRU_LEN *  CRU_WID;    % 0.37 [um/um^2] * ( 1.84 [um] * 0.9 [um] )[/CRU] -> [um/CRU]
    tubule_density_per_cru_2d_std   = EXP_DEN_STD * CRU_LEN *  CRU_WID;
    
%     id_file     = 1;
%     while id_file <= TUBULE_POPULATION_NUM
%         tubule_density_per_cru_2d  = normrnd(tubule_density_per_cru_2d_mean, tubule_density_per_cru_2d_std);    % [um/CRU]
%         if (tubule_density_per_cru_2d > 0)
%             tubule_length_set_2d            = tubule_density_per_cru_2d * inner_len * inner_wid * inner_dep;    % [um / CRU] -> [um / Cell]
%             tubule_length_set(id_file, 1)   = tubule_length_set_2d * (EXP_AT_TT + 1 + inner_dep/inner_wid) / (EXP_AT_TT + 1); 	% 2D -> 3D: add z-direction tubule
%             tubule_length_set(id_file, 1)   = tubule_length_set(id_file, 1) * (EXP_BRANCH_MEAN + EXP_BRANCH_MEAN / (CRU_LEN + CRU_WID) * (EXP_AT_DIAMETER_MEAN + EXP_TT_DIAMETER_MEAN) * 2 * 3) / EXP_BRANCH_MEAN; 
%             id_file = id_file + 1;
%         end
%     end
%     tubule_length_set   = sort(tubule_length_set(:,1),'descend');
    
    % the default tubular density data set
    tubule_length_set   = [4892.94946905416; 4780.00588094376; 4711.37501721348; 4456.66091301611; 4146.45709464259; ...
                       4079.12231432912; 3964.67475385481; 3848.97564392040; 3528.06934892807; 3338.93734555472; ...
                       2761.10340929164; 2500.50731880106; 2337.07547702255; 2253.19779043250; 1809.45671138218; ...
                       1341.16593277257; 1167.05635013773; 973.329460211456; 225.615450293407; 0];

    exp_converted_branch_mean   = EXP_BRANCH_MEAN * (EXP_AT_TT + 1 + inner_dep/inner_wid) / (EXP_AT_TT + 1);
    unit_branch_length_tmp      = CRU_LEN * EXP_AT_TT + CRU_WID * (1 + inner_dep/inner_wid);
    tubule_length_set           = tubule_length_set .* (exp_converted_branch_mean + exp_converted_branch_mean / unit_branch_length_tmp * (EXP_TT_DIAMETER_MEAN * (1 + inner_dep/inner_wid) * 2 + EXP_AT_DIAMETER_MEAN * EXP_AT_TT * 2 )) / exp_converted_branch_mean;                      
end