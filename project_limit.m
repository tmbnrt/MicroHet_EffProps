%%%%%%%%%%%%%%%%%%%%%%%%%% Purpose%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this file is to see how the curve of the analytical solution
% approaches the curve for the rigid inclusions when E_I/E_m -> INF, and
% the similarly for the case of void when E_I/E_m -> 0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% material properties
f = [10.0 5.0 1.0 0.25 0.1];         % E_I / E_M

E_m = 100;
ny_m = .3;
ny_i = .3;

kappa_m = E_m/(3*(1-2*ny_m));
mu_m    = E_m/(2*(1+ny_m));

ci = 0.000001:0.02:1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here you can select which approximation type you want to analyze by 
% setting code accordingly
% code = 1 ..... Dilute Distribution
% code = 2 ..... Mori-Tanaka
% code = 3 ..... Self-consistent
% code = 4 ..... Differential Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
code = 4

if 1
    h_mu = figure;
    set( 0, 'currentfigure', h_mu);
    h1_mu_DD = zeros( length(f));
    for j = 1:length(f)

        E_i = f(j) * E_m;
        kappa_i = E_i/(3*(1-2*ny_i));
        mu_i    = E_i/(2*(1+ny_i));
        % theoretical solutions by dilute distribution
        kappa_eff = zeros( length( ci), 1);
        mu_eff    = zeros( length( ci), 1);
        E_eff   = zeros( length( ci), 1);
        for i = 1:length(ci)
            switch code
                case 1
                  [ kappa_eff(i), mu_eff(i)] = dilute_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
                    'iso');
                case 2
                  [ kappa_eff(i), mu_eff(i)] = mtanaka_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
                    'iso');
                case 3
                  [ kappa_eff(i), mu_eff(i)] = sc_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
                    'iso');
                case 4
                  [ kappa_eff(i), mu_eff(i)] = diff_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
                    'iso');
            end
            E_eff(i) = 9*kappa_eff(i)*mu_eff(i)/(3*kappa_eff(i)+mu_eff(i))
        end

        
        hold on;
        h1_mu_eff(j) = plot( ci, mu_eff, '-');        
        hold off;
    end
    xlabel( 'volume fraction of inclusion')
    ylabel( '\mu_{eff}')
end




