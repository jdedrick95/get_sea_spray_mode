% get_sea_spray_mode(bsca_RGB, bsca_std, bsca_inst_std, PNSD, PNSD_D, D_op, PNSD_std, PNSD_N_std, PNSD_D_std)
%
%
% Retrieves lognormal fitting parameters of the sea spray mode using
% submicron size distributions and supermicron scattering.
%
% Inputs:
%
% !!! REQUIRES sea_spray_mie_table.mat TO BE IN DIRECTORY OF FUNCTION !!!
%
% bsca_RGB: supermicron scattering coefficients (Mm^-1) [450nm, 550nm,
% 700nm] (time dimension in rows)
%
% bsca_std: standard deviation of supermicron scattering coefficients 
% during temporal average (Mm^-1) 450nm, 550nm, 700nm] 
% (time dimension in rows)
%
% bsca_inst_std: instrument scattering error/uncertainty
%
% PNSD: submicron particle size distribution (cm^-3 µm^-1) (time dimension in columns)
%
% PNSD_D: submicron particle diameters (µm)
%
% D_op: overlap region for which to constrain Mie solutions; nominally D >
% ~ 0.4 µm)
%
% PNSD_std: standard deviation of particle size distribition during
% temporal average; resolved at each size bin (cm^-3 µm^-1)
% 
% PNSD_N_std: instrument concentration error/uncertinaty 
%
% PNSD_D_std: instrument sizing uncertainty
%
%
%
% Outputs:
%
% sea_spray_mode: sea spray mode fitting parameters [number, mean diameter,
% geometric standard deviation]
%
% sea_spray_mode_95: 95% confidence interval ranges of the sea spray mode
% fitting parameters. ROWS: number, mean diameter, geometric standard
% deviation. COLUMNS: lower 95th, upper 95th.
%
%
% test_coeff: probable Mie solutions that are tested against the size 
% distribution
% 
% RSS_fit: residual sum of squares of unique sea spray mode to measured
% size distribution
%
% chi2_fit: chi-square error of unique sea spray mode to measured size
% distribution
%
% D_mie: size distribution diameters (um)
%
% dlogDp_mie: log10 difference of the diameters
%
% fail_flag: flag value identifying reason for retrieval failure (0 =
% retrieval successful, 1 = scattering not available at all 3 wavelengths,
% 2 = no Mie scattering solutions below the error threshold, 3 = no Mie
% solutions that are within the joint probability 95th percentile that can
% be tested against the size distribution).
%
% retrieval_duration: time to complete the retrieval (minutes)
%
%
%
% Written by: Jeramy L. Dedrick
%             Scripps Institution of Oceanography
%             Last edited 16 January 2022
% 
% For description of methodology, see Dedrick et al. (2022), DOI: #



function [sea_spray_mode,...
          sea_spray_mode_95,...
          error_thresh,...
          low_error_idx,...
          test_coeff,...
          RSS_fit,...
          chi2_fit,...
          D_mie,...
          dlogDp_mie,...
          fail_flag,...
          retrieval_duration] = get_sea_spray_mode(bsca_RGB, ...
                                                   bsca_std, ...
                                                   bsca_inst_std, ...
                                                   PNSD, ...
                                                   PNSD_D, ...
                                                   D_op, ...
                                                   PNSD_std, ...
                                                   PNSD_N_std, ...
                                                   PNSD_D_std)





 
% loads look up table
load('sea_spray_mie_table.mat');


N_data     = size(bsca_RGB,1);              % size of data
D_mie      = D';                            % particle diameter (µm) (from Mie Table)
dlogDp_mie = nanmean(diff(log10(D_mie)));   % dlogDp (from Mie Table)
                      
 

%% selecting fitting region

fit_idx = [find(PNSD_D == D_op(1)), find(PNSD_D == D_op(end))];
fit_idx = fit_idx(1):fit_idx(2);
fit_D   = PNSD_D(fit_idx);

           


%% Running Retrieval 

for i = 1:N_data
              
    tic
    
    % FAIL CHECK: if scattering data not available at all wavelengths              
    fail_check_1 = sum(isnan(bsca_RGB(i,:)));


    if fail_check_1 > 0 % FAIL CODE, at least one wavelength missing
        
        fail_flag(i)                = 1;
        
        sea_spray_mode(i,:)         = NaN(1,3);
        sea_spray_mode_95{i}        = NaN(3,3);
        error_thresh(i)             = NaN(1,1);
        low_error_idx{i}            = NaN(1,1);
        test_coeff{i}               = NaN(1,3);  
        RSS_fit(i)                  = NaN(1,1);
        chi2_fit(i)                 = NaN(1,1);
        retrieval_duration(i)       = toc / 60;
        
              
    else
        
        % scatter error between NEPH and MIE
        dbsca_B = abs(bsca_RGB(i,1) - b_sca(:,1));
        dbsca_G = abs(bsca_RGB(i,2) - b_sca(:,2));
        dbsca_R = abs(bsca_RGB(i,3) - b_sca(:,3));
        
        % total error of scattering
        dbsca_RGB = sqrt( (dbsca_B).^2 + ...
                          (dbsca_G).^2 + ...
                          (dbsca_R).^2 );
        
        % error threshold  
        error_thresh(i) = sqrt( nansum( (bsca_std(i,:).^2 + bsca_inst_std.^2) ) );        
        
        % solutions below the error threshold
        low_error_idx{i}  = find(dbsca_RGB < error_thresh(i));
        
        % FAIL CHECK: if no solutions are below error threshold (fail code)
        fail_check_2      = isempty(low_error_idx{i}); 
        
        
        
        if fail_check_2 == 1
            
            fail_flag(i)                = 2;
            
            sea_spray_mode(i,:)         = NaN(1,3);
            sea_spray_mode_95{i}        = NaN(3,3);
            test_coeff{i}               = NaN(1,3);  
            RSS_fit(i)                  = NaN(1,1);
            chi2_fit(i)                 = NaN(1,1);
            retrieval_duration(i)       = toc / 60;            
            
            
        else
            
          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Reducing Sample Space to Most Probable Solutions %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % gets fitting parameters from solutions below error threshold
        coeff_get   = [coeff(low_error_idx{i},1),...
                       coeff(low_error_idx{i},2),...
                       coeff(low_error_idx{i},3)];

        % size of the reduced solution space (from error threshold)           
        N_coeff_get = size(coeff_get,1); 
     
        
        % ranges of retrieved fitting parameters
        % Number
        coeff_N_space    = [min(coeff_get(:,1)), ...                  % minimum
                            mean(diff(unique(coeff_get(:,1)))), ...   % interval/step
                            max(coeff_get(:,1))];                     % maximum

        % Geometeric Mean Diameter             
        coeff_Dg_space   = [min(coeff_get(:,2)), ...                  % minimum
                            mean(diff(unique(coeff_get(:,2)))), ...   % interval/step
                            max(coeff_get(:,2))];                     % maximum            

        % Geometric Standard Deviation             
        coeff_sigg_space = [min(coeff_get(:,3)), ...                  % minimum
                            mean(diff(unique(coeff_get(:,3)))), ...   % interval/step
                            max(coeff_get(:,3))];                     % maximum          
        
      
        
%%%%%%%%%%%%%%%%%%%% Number | Width joint probablity %%%%%%%%%%%%%%%%%%%%%%

        % computes joint probability
        [joint_N_sig_count_r, ...
         joint_N_sig_centers_r] = hist3([coeff_get(:,1), coeff_get(:,3)], ...
                                         'ctrs', ...
                                         {coeff_N_space(1):coeff_N_space(2):coeff_N_space(3) ...
                                          coeff_sigg_space(1):coeff_sigg_space(2):coeff_sigg_space(3)});
       
                                                    
        % joint probability (normalized)                                           
        joint_N_sig_prob_r = joint_N_sig_count_r ./ N_coeff_get; 
        joint_N_sig_prob_r(joint_N_sig_prob_r == 0) = NaN; % NaNs zero probability
        
        
        % grabs ranges of parameters
        N_sig_N   = joint_N_sig_centers_r{1,1}; % Y (rows)
        N_sig_sig = joint_N_sig_centers_r{1,2}; % X (columns)

        % computes joint probability 95th percentile 
        N_sig_prctile   = prctile(joint_N_sig_prob_r, 95, 'all');

        % finds where joint probability is >= 95th percentile
        [N_sig_prctile_r, N_sig_prctile_c] = find(joint_N_sig_prob_r >= N_sig_prctile);

        % creates grid of parameters to grab values above probability
        % percentile threshold
        [N_sig_x, N_sig_y] =  meshgrid(N_sig_sig, N_sig_N);


        % find indices of low error solutions that have high probability
        % combination 
        test_coeff_N_sig = [];
        for idx_N_sig = 1:length(N_sig_prctile_r)

            % parameter values from grid
            N_sig_mapx(idx_N_sig) = N_sig_x(N_sig_prctile_r(idx_N_sig),N_sig_prctile_c(idx_N_sig));
            N_sig_mapy(idx_N_sig) = N_sig_y(N_sig_prctile_r(idx_N_sig),N_sig_prctile_c(idx_N_sig));

            % low error solution indices
            test_coeff_N_sig_r = find(coeff_get(:,1) == N_sig_mapy(idx_N_sig) & ...
                                      coeff_get(:,3) == N_sig_mapx(idx_N_sig));

            % concatenated indices                    
            test_coeff_N_sig = [test_coeff_N_sig; test_coeff_N_sig_r];


        end         
        
        % fitting parameter solutions to test against size distribution
        test_coeff_r  = coeff_get(test_coeff_N_sig,:);
        test_coeff{i} = test_coeff_r;
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Fails code if there are no testable solutions
        if isempty(test_coeff_r)
            
            fail_flag(i) = 3;
            
            test_coeff{i}               = NaN;
            sea_spray_mode(i,:)         = NaN(1,3);
            sea_spray_mode_95{i}        = NaN(3,3);
            test_coeff{i}               = NaN(1,3); 
            RSS_fit(i)                  = NaN(1,1);
            chi2_fit(i)                 = NaN(1,1);
            retrieval_duration(i)       = toc / 60;
            
        else
        
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%% SIZE DISTRIBUTION PERTURBATION %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % test group size
        N_group(i)       = size(test_coeff_r,1);

        % number of times to perturb size distribution
        N_perturb        = 100;
        
        % variability/uncertainty in PNSD
        PNSD_uncertainty = sqrt(PNSD_std(:,i).^2 + PNSD_N_std.^2);
        
        % chi squared critical value
        chi2_crit        = chi2inv(0.95, length(fit_idx) - 1);
        
        
        % perturbing size distribution
        for j = 1 : N_perturb + 1
            
            if j == 1 % unperturbed 
                
                PNSD_perturb(:,j)   = PNSD(:, i); 
                diam_perturb(:,j)   = PNSD_D; 
                PMSD_perturb(:,j)   = PNSD_perturb(:,j) .* ...
                                      (pi / 6) .* ...
                                      (1) .* ...
                                      (diam_perturb(:,j) .* 1e-6) .^3 .* ...
                                      (1e6 .* 1e6) .* 1e6; 
                
            else 
                
                % perturbing distribution
                PNSD_perturb(:,j)   = normrnd(PNSD(:,i), PNSD_uncertainty);
                
                % zero out values that are less than zero                          
                zero_out_idx                 = find(PNSD_perturb(:,j) < 0);
                PNSD_perturb(zero_out_idx,j) = NaN;                          
                zero_out_idx                 = [];                         
                                          
                % perturbing size     
                diam_perturb(:,j) = PNSD_D + (PNSD_D .* normrnd(0,PNSD_D_std,[1,1])); 
             
                % converting PNSD to PMSD                          
                PMSD_perturb(:,j) = PNSD_perturb(:,j) .* ...
                                    (pi / 6) .* ...
                                    (1) .* ...
                                    (diam_perturb(:,j) .* 1e-6) .^3 .* ...
                                    (1e6 .* 1e6) .* 1e6;                          
                                          
                                          
            end % IF: perturbing condition (1 = obs/mean, >1 = perturb)
            
            
            for k = 1:N_group(i)                         
                
                % fitting parameters from low error cluster
                N_MIE    = test_coeff_r(k,1); % number concentration
                Dg_MIE   = test_coeff_r(k,2); % geometric mean diameter
                sig_MIE  = test_coeff_r(k,3); % geometric standard deviation
                
                % perturbed diameter
                diam_MIE = diam_perturb(:,j);
                
                
                % derived PNSD
                PNSD_MIE_deriv(:,k,j) = N_MIE ./ ...
                                        (sqrt(2*pi) .* log10(sig_MIE)) .* ...
                                        exp( - ((log10(diam_MIE) - log10(Dg_MIE)).^2) ./ ...
                                        (2 .* log10(sig_MIE).^2));                     
                
                % derived PMSD                    
                PMSD_MIE_deriv(:,k,j) = PNSD_MIE_deriv(:,k,j) .* ...
                                        (pi / 6) .* ...
                                        (1) .* ...
                                        (diam_MIE .* 1e-6) .^3 .* ...
                                        (1e6 .* 1e6) .* 1e6;                      
        
                                    
                % residual error of fit
                RSS_test_group(:,k,j)     = nansum( (PMSD_perturb(fit_idx,j) - PMSD_MIE_deriv(fit_idx,k,j)).^2 );

                % chi-square error
                chi2_test_group(:,k,j)    = nansum(((PMSD_MIE_deriv(fit_idx,k,j) - PMSD_perturb(fit_idx,j)).^2) ./ PMSD_perturb(fit_idx,j));
        
                % removes fits with a chi-square error above the critical
                % threshold
                if chi2_test_group(:,k,j) > chi2_crit
                    
                    chi2_test_group(:,k,j) = NaN;
                    
                else 
                    
                    chi2_test_group(:,k,j) = chi2_test_group(:,k,j);
                    
                end
                
                
            end % FOR, k: test group 
            
            
            % finds minimum RSS for each perturbed PMSD
            RSS_test_group_min(:,j)   = nanmin(RSS_test_group(:,:,j));
            
            % finds minimum chi-square for each perturbed PMSD
            chi2_test_group_min(:,j)  = nanmin(chi2_test_group(:,:,j));
            
            % index of minimum residual for each perturbation
            RSS_min_idx(:,j)       = find(RSS_test_group(:,:,j) == RSS_test_group_min(:,j),1);
            
            if isnan(chi2_test_group_min(:,j))
                
                chi2_min_idx(:,j) = 1;
                
            else
            
            % index of minimum chi-square for each perturbation
            chi2_min_idx(:,j)      = find(chi2_test_group(:,:,j) == chi2_test_group_min(:,j),1);
            
            end
            
            
            % fitting coefficients of minimum residual solution for each
            % perturbation
            coeff_group_min(j,:) = [test_coeff_r(RSS_min_idx(:,j),1), ...
                                    test_coeff_r(RSS_min_idx(:,j),2), ...
                                    test_coeff_r(RSS_min_idx(:,j),3)];
                                
           % fitting coefficients of minimum chi-square solution for each
           % perturbation
            coeff_group_min_chi2(j,:) = [test_coeff_r(chi2_min_idx(:,j),1), ...
                                         test_coeff_r(chi2_min_idx(:,j),2), ...
                                         test_coeff_r(chi2_min_idx(:,j),3)];                     
                                
        end % FOR, j: perturbing PNSD and diameter
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        


       % z value for 95% confidence interval
       t_stat95 = tinv([0.025 0.975], size(coeff_group_min,1));

       % 95% confidence interval of fitting parameters
       coeff_95_r = [nanmean(coeff_group_min(:,1)) + (t_stat95 .* (nanstd(coeff_group_min(:,1)) ./ sqrt(size(coeff_group_min,1) - 1))); ...
                     nanmean(coeff_group_min(:,2)) + (t_stat95 .* (nanstd(coeff_group_min(:,2)) ./ sqrt(size(coeff_group_min,1) - 1))); ...
                     nanmean(coeff_group_min(:,3)) + (t_stat95 .* (nanstd(coeff_group_min(:,3)) ./ sqrt(size(coeff_group_min,1) - 1)))];

       % 95% confidence interval of fitting parameters (chi-square)
       coeff_95_chi2_r = [nanmean(coeff_group_min_chi2(:,1)) + (t_stat95 .* (nanstd(coeff_group_min_chi2(:,1)) ./ sqrt(size(coeff_group_min_chi2,1) - 1))); ...
                          nanmean(coeff_group_min_chi2(:,2)) + (t_stat95 .* (nanstd(coeff_group_min_chi2(:,2)) ./ sqrt(size(coeff_group_min_chi2,1) - 1))); ...
                          nanmean(coeff_group_min_chi2(:,3)) + (t_stat95 .* (nanstd(coeff_group_min_chi2(:,3)) ./ sqrt(size(coeff_group_min_chi2,1) - 1)))];
                 
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%% Removing Solutions not within 95% CI %%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% floor and ciel functions used to account for resolution of parameters
       
       % number concentration solutions outside of probable range
       N95L_remove = find(coeff_group_min(:,1) < floor(coeff_95_r(1,1)) - coeff_N_space(2) / 2);
       N95U_remove = find(coeff_group_min(:,1) > ceil(coeff_95_r(1,2)) + coeff_N_space(2) / 2);
       Nremove     = unique([N95L_remove; N95U_remove]);
       
       % mean diameter solutions outside of probable range
       D95L_remove = find(coeff_group_min(:,2) < (floor(coeff_95_r(2,1) * 100) / 100) - coeff_Dg_space(2) / 2);
       D95U_remove = find(coeff_group_min(:,2) > (ceil(coeff_95_r(2,2) * 100) / 100) + coeff_Dg_space(2) / 2);
       Dremove     = unique([D95L_remove; D95U_remove]);
       
       % geometric standard devation solutions outside of probable range
       S95L_remove = find(coeff_group_min(:,3) < (floor(coeff_95_r(3,1) * 10) / 10) - coeff_sigg_space(2) / 2); 
       S95U_remove = find(coeff_group_min(:,3) > (ceil(coeff_95_r(3,2) * 10) / 10) + coeff_sigg_space(2) / 2);
       Sremove     = unique([S95L_remove; S95U_remove]);
              
       % selecting only the most probable solutions
       coeff_prob            = coeff_group_min;
       coeff_prob(Nremove,1) = NaN;
       coeff_prob(Dremove,2) = NaN;
       coeff_prob(Sremove,3) = NaN;
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
       % unique sea spray mode solution [number, mean diameter, width]
       sea_spray_mode(i,:) = nanmean(coeff_prob, 1);
              
       % 95% confidence interval of unique solution fitting parameters
       sea_spray_mode_95{i} = coeff_95_r;     

       
       % unique sea spray mode size distribution
       sea_spray_mode_PNSD_fit = sea_spray_mode(i,1) ./ ...
                                 (sqrt(2*pi) .* log10(sea_spray_mode(i,3))) .* ...
                                 exp( - ((log10(fit_D) - log10(sea_spray_mode(i,2))).^2) ./ ...
                                 (2 .* log10(sea_spray_mode(i,3)).^2));    
                                    
       sea_spray_mode_PMSD_fit = sea_spray_mode_PNSD_fit .* ...
                                 (pi / 6) .* ...
                                 (1) .* ...
                                 (fit_D .* 1e-6) .^3 .* ...
                                 (1e6 .* 1e6) .* 1e6;  
                             
        % measured PMSD
        PMSD_fit = PNSD(fit_idx, i) .* ...
                   (pi / 6) .* ...
                   (1) .* ...
                   (fit_D .* 1e-6) .^3 .* ...
                   (1e6 .* 1e6) .* 1e6;  
                             
        % residual error of unique sea spray mode
        RSS_fit(i)    = nansum( (PMSD_fit - sea_spray_mode_PMSD_fit).^2 );

        % chi-square error of unique sea spray mode
        chi2_fit(i)   = nansum(((PMSD_fit - sea_spray_mode_PMSD_fit).^2) ./ PMSD_fit);


       
       fail_flag(i) = 0;
       
       
       
       % clearing variables for future use
          clear PNSD_perturb ...        
          PNSD_std_perturb ...
          diam_perturb ...
          PMSD_perturb ...
          diam_MIE ...
          PNSD_MIE_deriv ...
          PMSD_MIE_deriv ...
          RSS_test_group ...
          RSS_test_group_min ...
          RSS_min_idx ...
          chi2_test_group ...
          chi2_test_group_min ...
          chi2_min_idx ...
          coeff_group_min ...
          coeff_group_min_chi2 ...
          t_stat95 ...
          coeff_95_r ...
          coeff_95_chi2_r ...
          MIE_deriv_PNSD_1 ...
          MIE_deriv_PMSD_1 ...
          PNSD_uncertainty ...
          coeff_prob ...
          Nunique ...
          Dunique ...
          Sunique ...
          sea_spray_mode_PNSD_fit ...
          sea_spray_mode_PMSD_fit ...
          PMSD_fit ...
          dbsca
      
       

        end % Fail Check 3
        
        
        
        end % Fail Check 2
        
  
        
    end % Fail Check 1
    
    
    
    retrieval_duration(i) = toc / 60;
    
    
    disp(strcat('Completed Retrieval:', num2str(i),'/', num2str(N_data)))
    
end % loop through retrievals
        
        
        
end % function
