
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This m.file is employed for running the simulations of the astro-field model to produce
% the data and figures corresponding to the results presented in Figures
% 3-7 of the article
% 
% [1] An astro-neural-field model with application to cortical spreading
% depolarization, E. Baspinar, D. Avitabile, C. Nouveau, M. Desroches, F. Campillo, M.
% Mantegazza.
%
% This code was designed and edited by
% Daniele Avitabile, VU Amsterdan, Netherlands
% d.avitabile@vu.nl
% Emre Baspinar, Inria, MathNeuro Team, France
% emre.baspinar@inria.fr 
%
% It was contributed by F. Campillo.
%
% Contact: emre.baspinar@inria.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Run the m.files by index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     
%%%     index = 1 --> ISO_or_GBZ.m            (Figure 3a in [1])
%%%     index = 2 --> ChR2.m                  (Figure 3b in [1]) 
%%%     index = 3 --> GBZ_plus_ChR2.m         (Figure 3c in [1])
%%%     index = 4 --> sweep_gamma_aa.m        (Figure 4  in [1])
%%%     index = 5 --> sweep_gamma_a.m         (Figure 5  in [1])
%%%     index = 6 --> beyond_propagation.m    (Figure 6  in [1])
%%%     index = 7 --> astrocyte_stops.m       (Figure 7  in [1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reset environment
clear all; close all; clc;

% ---- SET INDEX HERE ----
index = 7;   % choose 1, 2, 3, 4, 5, 6 or 7, nothing else.
    
indices = [1, 2, 3, 4, 5, 6, 7];

% ---- CHECK VALIDITY ----
if ~ismember(index, indices)

    error('Invalid index %d. Available indices are: %s', index, mat2str(indices));

elseif index == 1  

    script_name = 'ISO_or_GBZ';

    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);

elseif index == 2

    script_name = 'ChR2';
    
    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);

elseif index == 3

    script_name = 'GBZ_plus_ChR2';
    
    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);

elseif index == 4

    script_name = 'sweep_gamma_aa';
    
    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);

elseif index == 5

    script_name = 'sweep_gamma_a';
    
    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);

elseif index == 6

    script_name = 'beyond_propagation';
    
    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);
    
elseif index == 7

    script_name = 'astrocyte_stops';
    
    % ---- RUN SCRIPT ----
    fprintf('Running %s.m...\n', script_name);
    run(script_name);

end
