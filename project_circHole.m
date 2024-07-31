%%%%%%%%%%%%%%%%%%%%%%%%%% Purpose%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Considering the 2D-example of a plate in plane stress with circular
% holes, see an interesting characteristic of the self-consistent method
% that can predict the total loss of the macroscopic stiffness below ci->1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clc

% material properties of the plate
E = 100;
nu = .3;
mu = E/(2*(1+nu));

% area fraction of the holes
ci = 0:0.02:1;

% initialize Elastic Young's modulus and shear modulus
% dilute distribution
E_DD = zeros( length( ci), 1);
mu_DD = zeros( length( ci), 1);

% Mori-Tanaka distribution
E_MT    = zeros( length( ci), 1);
mu_MT   = zeros( length( ci), 1);

% self-consistent method
E_SC    = zeros( length( ci), 1);
mu_SC  	= zeros( length( ci), 1);

% differential scheme
E_DF    = zeros( length( ci), 1);
mu_DF  	= zeros( length( ci), 1);


% loop over all area fractions of the holes
for i = 1:length(ci)

    % dilute distribution
    E_DD(i)     = E/(1+3*ci(i));
    mu_DD(i)    = E/(2*(1+nu+4*ci(i)));

    % Mori-Tanaka distribution
    E_MT(i)     = E*(1-ci(i))/(1+2*ci(i));
    mu_MT(i) 	= mu*((1-ci(i))*(1+nu))/(1+nu+ci(i)*(3-nu));

    % self-consistent method
    E_SC(i)     = E*(1-3*ci(i));
    mu_SC(i)  	= E*(1-3*ci(i))/( 2*(1+ci(i)+nu*(1-3*ci(i))) );

  % differential scheme
    E_DF(i)     = E*(1-ci(i))^3;
    mu_DF(i)    = mu*( 3*(1+nu)*(1-ci(i))^3 ) / (4+(3*nu-1)*(1-ci(i))^3);

end

% plot E
figure('Name', 'E*')
title('E^*', 'Fontsize', 16)
xlabel('area fraction of circular holes', 'Fontsize', 16)
ylabel('Young''s modulus [MPa]', 'Fontsize', 16)
hold on
plot( ci, E_DD, 'o-r');
plot( ci, E_MT, '.-g');
plot( ci, E_SC, 'x-b');
plot( ci, E_DF, '.-k');
hold off
legend('DD','MT','SC','DS');
set(gca, 'YLim', [0 E])
box on

% plot mu
figure('Name', 'mu^*')
title('\mu^*', 'Fontsize', 16)
xlabel('area fraction of circular holes', 'Fontsize', 16)
ylabel('Shear modulus [MPa]', 'Fontsize', 16)
hold on
plot( ci, mu_DD, 'o-r');
plot( ci, mu_MT, '.-g');
plot( ci, mu_SC, 'x-b');
plot( ci, mu_DF, '.-k');
hold off
legend('DD','MT','SC','DS');
set(gca, 'YLim', [0 mu])
box on



