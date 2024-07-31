%%%%%%%%%%%%%%%%%%%%%%%%%% Purpose%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this file shows you how to compute the theoretical solutions for 
% all the analytical methods and provides graphical illustrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clc

% material properties
Ee_m  = 47;
Ee_i  = 470;
nuu_m = 0.29;
nuu_i = 0.3;

kappa_m =    (Ee_m)/(3*(1-2*(nuu_m)));
mu_m    =    (Ee_m)/(2*(1+(nuu_m)));
kappa_i =    (Ee_i)/(3*(1-2*(nuu_i)));
mu_i    =    (Ee_i)/(2*(1+(nuu_i)));

% volume fraction of inclusion
ci = 0:0.02:1;

% Voigt bound
kappa_Voigt  = ci * kappa_i  + ( 1 - ci) * kappa_m;
mu_Voigt     = ci * mu_i + ( 1 - ci) * mu_m;

% Reuss bound
kappa_Reuss  = ci/kappa_i  + ( 1 - ci)/kappa_m;
kappa_Reuss  = 1./kappa_Reuss;

mu_Reuss     = ci/mu_i + ( 1 - ci)/mu_m;
mu_Reuss     = 1./mu_Reuss;

% plot Voigt- and Reuss bound
h_kappa = figure;
h1_kappa_V = plot(ci, kappa_Voigt, 'g');
hold on
h1_kappa_R = plot(ci, kappa_Reuss, 'b');
hold off
title('\kappa^*', 'Fontsize', 16)
xlabel('volume fraction of inclusion', 'Fontsize', 16)
ylabel('bulk modulus [MPa]', 'Fontsize', 16)
legend('Voigt','Reuss');

h_mu = figure;
h1_mu_V = plot(ci, mu_Voigt, 'g');
hold on
h1_mu_R = plot(ci, mu_Reuss, 'b');
hold off
title('\mu^*', 'Fontsize', 16)
xlabel('volume fraction of inclusion', 'Fontsize', 16)
ylabel('shear modulus [MPa]', 'Fontsize', 16)
legend('Voigt','Reuss');

% Dilute-Distribution Method
if 1
    kappa_DD = zeros( length( ci), 1);
    mu_DD    = zeros( length( ci), 1);
    for i = 1:length(ci)
        [ kappa_DD(i), mu_DD(i)] = dilute_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
            'iso');
    end

    % plot
    set( 0, 'currentfigure', h_kappa);
    hold on;
    h1_kappa_DD = plot( ci, kappa_DD, 'o-r');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution');
    hold off;

    set( 0, 'currentfigure', h_mu);
    hold on;
    h1_mu_DD = plot( ci, mu_DD, 'o-r');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution', 'Location', 'NorthWest');
    hold off;
end

 % Mori-Tanaka Approach
if 1
    kappa_MT = zeros( length( ci), 1);
    mu_MT    = zeros( length( ci), 1);
    for i = 1:length(ci)
        [ kappa_MT(i), mu_MT(i)] = mtanaka_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
            'iso');
    end

    % plot
    set( 0, 'currentfigure', h_kappa);
    hold on;
    h1_kappa_MT = plot( ci, kappa_MT, '.-r');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution', 'Mori-Tanaka', 'Location', 'NorthWest');
    hold off;

    set( 0, 'currentfigure', h_mu);
    hold on;
    h1_mu_MT = plot( ci, mu_MT, '.-r');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution', 'Mori-Tanaka', 'Location', 'NorthWest');
    hold off;
end

% Self-Consistent Method
if 1
    kappa_SK = zeros( length( ci), 1);
    mu_SK    = zeros( length( ci), 1);

    for i = 1:length(ci)
        [ kappa_SK(i), mu_SK(i)] = sc_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
            'iso');
    end

    % plot
    set( 0, 'currentfigure', h_kappa);
    hold on;
    h1_kappa_SK = plot( ci, kappa_SK, 'x-m');
    legend( 'Voigt', 'Reuss', 'DD', 'MT', 'SC', 'Location', 'NorthWest');
    hold off;

    set( 0, 'currentfigure', h_mu);
    hold on;
    h1_mu_SK = plot( ci, mu_SK, 'x-m');
    legend( 'Voigt', 'Reuss', 'DD', 'MT', 'SC', 'Location', 'NorthWest');
    hold off;
end

% Differential Scheme
if 1
    kappa_DF = zeros( length( ci), 1);
    mu_DF    = zeros( length( ci), 1);

    % c_i = 0 and c_i = 1 can't be computed by the Matlab server
    kappa_DF(1) = kappa_m; kappa_DF(end) = kappa_i;
    mu_DF(1)    = mu_m;    mu_DF(end)    = mu_i;

    for i = 2:length(ci)-1
        [ kappa_DF(i), mu_DF(i)] = diff_analy( ci(i), kappa_m, mu_m, kappa_i, mu_i, ...
            'iso');
    end

    % plot
    set( 0, 'currentfigure', h_kappa);
    hold on;
    h1_kappa_DF = plot( ci, kappa_DF, '.-y');
    legend( 'Voigt', 'Reuss', 'DD', 'MT', 'SC', 'DS', 'Location', 'NorthWest');
    hold off;

    set( 0, 'currentfigure', h_mu);
    hold on;
    h1_mu_DF = plot( ci, mu_DF, '.-y');
    legend( 'Voigt', 'Reuss', 'DD', 'MT', 'SC', 'DS', 'Location', 'NorthWest');
    hold off;
end


% Hashin-Shtrikman Bounds
if 0
    kappa_HS1  = ( 1/(kappa_i-kappa_m) + 3*(1-ci)/(3*kappa_m+4*mu_m) );
    kappa_HS1  = kappa_m + ci*1./(kappa_HS1);
    
    kappa_HS2  = ( 1/(kappa_m-kappa_i) + 3*ci/(3*kappa_i+4*mu_i) );
    kappa_HS2  = kappa_i + (1-ci)*1./(kappa_HS2);
    
    mu_HS1 = ( 1/(mu_i-mu_m) + 6*(1-ci)*(kappa_m+2*mu_m)/(5*mu_m*(3*kappa_m+4*mu_m)) );
    mu_HS1 = mu_m + ci*1./(mu_HS1);
    
    mu_HS2 = ( 1/(mu_m-mu_i) + 6*ci*(kappa_i+2*mu_i)/(5*mu_i*(3*kappa_i+4*mu_i)) );
    mu_HS2 = mu_i + (1-ci)*1./(mu_HS2);

    % plot
    set( 0, 'currentfigure', h_kappa); 
    hold on
    h1_kappa_HS1 = plot(ci, kappa_HS1, '--k');         
    h1_kappa_HS2 = plot(ci, kappa_HS2, ':k');
    hold off
    legend( 'Voigt', 'Reuss', 'DD', 'MT', 'SC', 'DS', 'HS-', 'HS+', 'Location', 'NorthWest');

    set( 0, 'currentfigure', h_mu); 
    hold on
    h1_mu_HS1 = plot(ci, mu_HS1, '--k');
    h1_mu_HS2 = plot(ci, mu_HS2, ':k');
    hold off
    legend( 'Voigt', 'Reuss', 'DD', 'MT', 'SC', 'DS', 'HS-', 'HS+', 'Location', 'NorthWest');

end


