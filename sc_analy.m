function [ kappa_eff, mu_eff] = sc_analy(c_i, kappa_m, mu_m, kappa_i, mu_i, str)
% c_i           = volume fraction of inclusion
% kappa_m, mu_m = kappa, mu of matrix phase
% kappa_i, mu_i = kappa, mu of inclusion phase

% tolerance
tol = 1e-6;
switch str
    case 'iso'
        % initial value
        kappa = [(1-c_i)*kappa_m + c_i*kappa_i;
                 (1-c_i)*mu_m    + c_i*mu_i];      %  !!!! error?!
        r     = [ 1; 1];
        
        % newton iteration
        while norm( r, 2) > tol
            % compute residual vector using function "get_residual"
            r = get_residual(c_i,[kappa(1) kappa_m kappa_i],[kappa(2) mu_m mu_i]);
            
            % compute tangent matrix using function "get_tangent"
            A = get_tangent(c_i,[kappa(1) kappa_m kappa_i],[kappa(2) mu_m mu_i]);
            
            % solve linear equation
            B = inv(A);
            dkappa = B*(-r);
            
            % update
            kappa  = kappa + dkappa ;
            
            
        end
        
        % effective compression modulus
        kappa_eff = kappa(1);
        
        % effective shear modulus
        mu_eff    = kappa(2);
        
    case 'plane strain'
end

function r = get_residual(c, k, u)
% input variables
% k  = [kappa_eff from last iter. step, kappa of matrix phase, kappa of inclusion phase]
% u  = [mu_eff from last iter. step,    mu of matrix phase,    mu of inclusion phase]

% volume fraction of the matrix phase
c2 = 1 - c;

% residual vector = [r1; r2]
r = [ 4*k(1)*u(1) + 3*(c2*k(3) + c*k(2))*k(1) - 4*(c2*k(2) + c*k(3))*u(1) - 3*k(3)*k(2); ...
      9*k(1)*u(1)^2 + 8*u(1)^3 + 3*k(1)*u(1)*((2-5*c)*u(3) + (2-5*c2)*u(2)) + ...
         + 4*u(1)^2*((3-5*c2)*u(2) + (3-5*c)*u(3)) - 6*k(1)*u(3)*u(2) ...
         - 12*u(1)*u(3)*u(2) ];


function A = get_tangent(c, k, u)
% input variables
% k  = [kappa_eff from last iter. step, kappa of matrix phase, kappa of inclusion phase]
% u  = [mu_eff from last iter. step,    mu of matrix phase,    mu of inclusion phase]

% volume fraction of the matrix phase
c2 = 1 - c;

% tangent matrix = [dr1/dk*  dr1/du*]
%                  [dr2/dk*  dr2/du*]
A = [ 4*u(1) + 3*(c2*k(3) + c*k(2)), ...
      4*k(1) - 4*(c2*k(2) + c*k(3)); ...
      9*u(1)^2 + 3*u(1)*((2-5*c)*u(3) + (2-5*c2)*u(2)) - 6*u(3)*u(2),...
      18*k(1)*u(1) + 24*u(1)^2 + 3*k(1)*((2-5*c)*u(3) + (2-5*c2)*u(2)) + 8*u(1)*((3-5*c2)*u(2) + (3-5*c)*u(3)) - 12*u(3)*u(2)];

