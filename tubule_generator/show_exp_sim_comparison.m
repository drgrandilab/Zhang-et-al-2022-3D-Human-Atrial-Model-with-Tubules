function show_exp_sim_comparison(CELL_LEN, CELL_WID, CELL_DEP, TUBULE_POPULATION_NUM)
    close all
    
    data_exp_rabbit = load('./EXP_RABBIT.mat');
    data_exp_human  = load('./EXP_HUMAN_DENSITY_BRANCH.mat');
    
    data_sim_central_plane_long_branch  = load('./validation_based_on_tubule_structure_single_plane_involve_z_branch.mat');
        
    population_size = TUBULE_POPULATION_NUM;
    
    x1_rand     = ones(population_size, 1) .* (1+(rand(population_size, 1)-0.5)/5);
    x2_rand     = ones(population_size, 1) .* (2+(rand(population_size, 1)-0.5)/5);

    x2  = 2.*ones(population_size, 1);

    x1  = 1.*ones(9, 1);
    figure(18)
    x   = [x1; x2];
    y   = [data_exp_human.EXP_HUMAN_DENSITY; data_sim_central_plane_long_branch.validate_tubule_density_inner_involve_z_dot(1:population_size)];
    f1  = boxplot(y, x, 'whisker', 0.7193, 'symbol', '', 'Color', 'k');
    set(f1, {'linew'}, {2})
    set(gca,'xcolor','none')
    hold on
    scatter(x2_rand, data_sim_central_plane_long_branch.validate_tubule_density_inner_involve_z_dot(1:population_size, 1), 'b', 'filled');
    hold off
    title('validation: exp v.s. sim')
    ylabel('Tubule Density (um/um^2)')
    set(gca, 'linew', 2, 'Fontsize', 20, 'TickDir', 'out')
    ylim([0 1])
    saveas(gcf,'Tubular_Density.pdf')
    saveas(gcf,'Tubular_Density.eps')
    
    figure(19)
    y3  = [data_sim_central_plane_long_branch.validate_tubule_at_branch(data_sim_central_plane_long_branch.validate_tubule_at_branch ~= 0); data_sim_central_plane_long_branch.validate_tubule_tt_branch(data_sim_central_plane_long_branch.validate_tubule_tt_branch ~= 0); data_sim_central_plane_long_branch.validate_tubule_z_tt_branch(data_sim_central_plane_long_branch.validate_tubule_z_tt_branch ~= 0)];
    x3  = 2.*ones(length(y3), 1);
    x   = [x1; x3];
    y   = [data_exp_human.EXP_HUMAN_BRANCH_LENGTH; y3];
    f1  = boxplot(y, x, 'whisker', 0.7193, 'symbol', '', 'Color', 'k');
    set(f1, {'linew'}, {2})
    set(gca,'xcolor','none')
    title('validation: exp v.s. sim')
    ylabel('Branch Length (um)')
    set(gca, 'linew', 2, 'Fontsize', 20, 'TickDir', 'out')
    ylim([0 6])
    saveas(gcf,'Branch_Length.pdf')
    saveas(gcf,'Branch_Length.eps')
    
    
    get_new_Cmem(CELL_LEN, CELL_WID, CELL_DEP, population_size)
    data_sim_Cmem   = load('get_new_Cmem.mat');
    
    figure(25)
    y3  = [data_sim_Cmem.Cmem_array];
    x   = 1.*ones(population_size, 1);
    y   = y3 .* 1e12;
    f1  = boxplot(y, x, 'whisker', 0.7193, 'symbol', '', 'Color', 'k');
    set(f1, {'linew'}, {2})
    set(gca,'xcolor','none')
    hold on
    scatter(x1_rand, y, 'b', 'filled');
    hold off
    title('validation: exp v.s. sim')
    ylabel('C_{mem} (pF)')
    set(gca, 'linew', 2, 'Fontsize', 20, 'TickDir', 'out')
    saveas(gcf,'Cmem.pdf')
    saveas(gcf,'Cmem.eps')
end