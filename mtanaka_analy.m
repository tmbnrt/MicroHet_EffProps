function [ kappa_eff, mu_eff] = mtanaka_analy( ci, kappa_m, mu_m, kappa_i, mu_i, str)
% obtaining the analytical effective material properties by Mori-Tanaka
% method
% ci:   volume fraction of the inclusion
% kappa_m:  kappa of the matrix material
% mu_m:     mu of the matrix material
% kappa_i:  kappa of the inclusion
% mu_i:     mu of the inclusion
switch str
    case 'iso'
        alpha = 3 * kappa_m / ( 3 * kappa_m + 4 * mu_m);
        beta  = 6 * ( kappa_m + 2 * mu_m) / 5 / ( 3 * kappa_m + 4 * mu_m);
        kappa_eff = kappa_m + ci * ( kappa_i - kappa_m) * kappa_m / ( kappa_m + alpha * ( 1 - ci) * ( kappa_i - kappa_m));
        mu_eff = mu_m + ci * ( mu_i - mu_m) * mu_m / ( mu_m + beta * ( 1 - ci) * ( mu_i - mu_m));
    case 'plane strain'
        E_m = 9 * kappa_m * mu_m / (3 * kappa_m + mu_m);
        E_i = 9 * kappa_i * mu_i / (3 * kappa_i + mu_i);
        ny_m = ( 3 * kappa_m - 2 * mu_m) / ( 6 * kappa_m + 2 * mu_m);
        ny_i = ( 3 * kappa_i - 2 * mu_i) / ( 6 * kappa_i + 2 * mu_i);

        C_i = get_C( E_i, ny_i);
        C_m = get_C( E_m, ny_m);
        S_m = get_S( ny_m);

        C_eff = C_m + ci * ( C_i - C_m) * inv( eye(3) + ( 1 - ci) * S_m * inv( C_m) * ( C_i - C_m));

        %         mu_eff = C_eff( 3, 3);
        lambda = C_eff( 1, 2);
        param_a = C_eff( 1, 1);
        ny_eff  = 1 / ( param_a /lambda + 1);
        E_eff  = param_a * ( 1 + ny_eff) *( 1 - 2 * ny_eff) / ( 1 - ny_eff);
        kappa_eff = E_eff / ( 3 * ( 1 - 2 * ny_eff));
        mu_eff    = E_eff / ( 2 * ( 1 + ny_eff));
        %         kappa_eff = lambda + 2 / 3 * mu_eff;
end

function C = get_C(E, ny)
C = E * ( 1 - ny) / ( 1 + ny) / ( 1 - 2 * ny) * ...
    [ 1, ny / ( 1 - ny), 0;
    ny / ( 1 - ny), 1, 0;
    0, 0, ( 1 - 2 * ny) / 2 / ( 1 - ny)];

function S = get_S( ny)

S = 1 / 8 / ( 1 - ny) * ...
    [ 5 - 4 * ny, 4 * ny - 1, 0;
    4 * ny - 1, 5 - 4 * ny, 0;
    0, 0, 3 - 4 * ny];