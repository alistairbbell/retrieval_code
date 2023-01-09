% function [quality_th,conv_th] = aura_mls_quality_parameters(species,vers)
%
% This function returns the quality and convergence threshold values and 
% the useful pressure range for different species retrieved by Aura MLS
% as given in the official Aura MLS Documentation.
% The parameters differ among the Level2 Retrieval versions!
%
% species = MLS Level2 species [string]
% vers = MLS Level2 Retrieval version (2.2 or 3.3)
%
% Example: [quality_th, conv_th, p_range] = aura_mls_quality_parameters('H2O',3.3);
%
% Dominik Scheiben, 2011-01-28
%
function [quality_th,conv_th,p_range] = aura_mls_quality_parameters(species,vers)

switch vers
    case 2.2
        switch upper(species)
            case 'H2O'
                quality_th = 0.9;
                conv_th = inf;
                p_range = [0.002 316];
            case 'O3'
                quality_th = 0.4;
                conv_th = 1.8;
                p_range = [0.022 215];
            case 'T'
                quality_th = 0.6;
                conv_th = 1.2;
                p_range = [0.001 316];
            case 'GPH'
                quality_th = 0.6;
                conv_th = 1.2;
                p_range = [0.001 316];
            case 'CLO'
                quality_th = 0.8;
                conv_th = 1.5;
                p_range = [1 100];
            case 'CO'
                quality_th = 0.2;
                conv_th = 1.8;
                p_range = [0.0046 215];
            case 'HNO3'
                quality_th = 0.4;
                conv_th = 1.8;
                p_range = [3.2 215];
            case 'BRO'
                quality_th = 1.2;
                conv_th = 1.5;
                p_range = [3.2 10];
            case 'HCL'
                quality_th = 1.0;
                conv_th = 1.5;
                p_range = [0.15 100];
            case 'HCN'
                quality_th = 0.2;
                conv_th = 2.0;
                p_range = [0.1 10];
            case 'HO2'
                quality_th = -inf;
                conv_th = 1.1;
                p_range = [0.032 22];
            case 'HOCL'
                quality_th = 1.4;
                conv_th = 1.5;
                p_range = [2.2 10];
            case 'N2O'
                quality_th = 0.5;
                conv_th = 1.55;
                p_range = [1 100];
            case 'OH'
                quality_th = -inf;
                conv_th = 1.1;
                p_range = [0.0032 32];
            case 'RHI'
                quality_th = 0.9;
                conv_th = inf;
                p_range = [0.002 316];
            case 'SO2'
                quality_th = 0.4;
                conv_th = 1.8;
                p_range = [10 215];
            case 'IWC'
                quality_th = -inf;
                conv_th = inf;
                p_range = [68 261];
            case 'IWP'
                quality_th = -inf;
                conv_th = inf;
                p_range = [0.001 1000];
            otherwise
                error('Unknown Aura MLS species!');
        end
    case 3.3
        switch upper(species)
            case 'H2O'
                quality_th = 1.3;
                conv_th = 2;
                p_range = [0.002 83];
            case 'O3'
                quality_th = 0.6;
                conv_th = 1.18;
                p_range = [0.02 100];
            case 'T'
                quality_th = 0.65;
                conv_th = 1.2;
                p_range = [0.001 261];
            case 'GPH'
                quality_th = 0.65;
                conv_th = 1.2;
                p_range = [0.001 261];
            case 'CLO'
                quality_th = 1.3;
                conv_th = 1.05;
                p_range = [1 147];
            case 'CO'
                quality_th = 0.2;
                conv_th = 1.4;
                p_range = [0.0046 100];
            case 'HNO3'
                quality_th = -inf;
                conv_th = inf;
                p_range = [1.5 215];
            case 'BRO'
                quality_th = 1.3;
                conv_th = 1.05;
                p_range = [3.2 10];
            case 'HCL'
                quality_th = 1.2;
                conv_th = 1.05;
                p_range = [0.32 100];
            case 'HCN'
                quality_th = 0.2;
                conv_th = 2.0;
                p_range = [0.1 10];
            case 'HO2'
                quality_th = -inf;
                conv_th = 1.1;
                p_range = [0.046 22];
            case 'HOCL'
                quality_th = 1.2;
                conv_th = 1.05;
                p_range = [2.2 10];
            case 'N2O'
                quality_th = 1.4;
                conv_th = 1.01;
                p_range = [0.46 100];
            case 'OH'
                quality_th = -inf;
                conv_th = 1.1;
                p_range = [0.0032 32];
            case 'RHI'
                quality_th = -inf;
                conv_th = inf;
                p_range = [0.002 316];
            case 'SO2'
                quality_th = 0.6;
                conv_th = 1.8;
                p_range = [10 215];
            case 'IWP'
                quality_th = -inf;
                conv_th = inf;
                p_range = [0.001 1000];
            case 'IWC'
                quality_th = -inf;
                conv_th = inf;
                p_range = [83 215];
            case 'CH3CL'
                quality_th = 1.3;
                conv_th = 1.05;
                p_range = [4.6 147];
            case 'CH3CN'
                quality_th = 1.4;
                conv_th = 1.05;
                p_range = [1 46];
            otherwise
                error('Unknown Aura MLS species!');
        end
    otherwise
        error('Unknown Aura MLS retrieval version!');
end