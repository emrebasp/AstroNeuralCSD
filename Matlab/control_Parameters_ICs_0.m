clear all; close all;



%% Parameters
    %% Model Parameters

    % Parameters of population transfer functions Se, Si, Sa    
    modelp.theta1_e = 1.7;              % theta1_e: initial potassium level of hyperactivation of the exc. population
    modelp.theta2_e = 1.8;              % theta2_e: initial potassium level of depolarization block of the exc. population
    modelp.theta1_i = 1.2;              % theta1_i: initial potassium level of hyperactivation of the inh. population
    modelp.theta2_i = 1.4;              % theta2_i: initial potassium level of depolarization block of the inh. population
    modelp.theta1_a = 1.0;              % theta1_a: initial potassium level of increased potassium clearance of the astro. population
    modelp.theta2_a = 1.2;              % theta2_a: initial potassium level of no potassium clearance  
    modelp.theta3_e  = 100;             % theta3_e: sharpness of the nonlinearity of S_e 
    modelp.theta3_i  = 10;              % theta3_i: sharpness of the nonlinearity of S_i 
    modelp.theta3_a  = 100;             % theta3_a: sharpness of the nonlinearity of S_a 
    modelp.theta4  = 0.3;               % theta4: threshold of the transfer functions S_e, S_i, S_a

    % Connection weights
    modelp.gamma_ei = 4;                % gamma_ei: connection weight from inh. to exc. neurons
    modelp.gamma_ie = 0.2;              % gamma_ie: connection weight from exc. to inh. neurons
    modelp.gamma_ee = 1.0;              % gamma_ee: connection weight from exc. to exc. neurons
    modelp.gamma_ii = 4;                % gamma_ii: connection weight from inh. to inh. neurons
    modelp.gamma_an = 0.5;              % gamma_an: connection weight from all neurons to astrocytes
    modelp.gamma_aa = 0.5;              % gamma_aa: connection weight from astrocytes to astrocytes
    
    % Parameters of potassium-dependent activation functions Fe, Fi, Fa
    modelp.alpha1_e  = 100;             % alpha1_e: sharpness of the nonlinerarity of Fe in ve equation 
    modelp.alpha1_i  = 10;              % alpha1_i: sharpness of the nonlinerarity of Fi in vi equation
    modelp.alpha1_a  = 100;             % alpha1_a: sharpness of the nonlinerarity of Fa in ra equation
    modelp.alpha2_e = 1.05;             % alpha2_e: potassium threshold of Fe function appearing in ve equation
    modelp.alpha2_i = 1.05;             % alpha2_i: potassium threshold of Fi function appearing in vi equation
    modelp.alpha2_a = 1.05;             % alpha2_a: potassium threshold of Fa function appearing in ra equation

    % Weights of potassium-dependent activation functions Fe, Fi, Fa
    modelp.gamma_eK = 1;         % weight for Fe
    modelp.gamma_iK = 1;         % weight for Fi
    modelp.gamma_aK = 1;         % weight for Fa

    % Parameters of potassium source function rho
    modelp.beta1_e = 10;                % beta1_e: sharpness of the nonlinearity of the rho function weighted by gamma_e in potassium equation
    modelp.beta2_e = 1.5;               % beta2_e: threshold value for the rho function weighted by gamma_e in potassium equation
    modelp.beta1_i = 10;                % beta1_i: sharpness of the nonlinearity of the rho function weighted by gamma_i in potassium equation
    modelp.beta2_i = 1.5*0.9975;        % beta2_i: threshold value for the rho function weighted by gamma_i in potassium equation
    modelp.beta1_a = 10;                % beta1_i: sharpness of the nonlinearity of the rho function weighted by gamma_a in potassium equation
    modelp.beta2_a = 1.5;               % beta2_a: threshold value for the rho function weighted by c3 in potassium equation

    % Weights of potassium source function rho
    modelp.gamma_e = 18.6;              % gamma_e: weight of contribution of exc. population spiking to extracellular potassium accumulation
    modelp.gamma_i = 18.6;              % gamma_i: weight of contribution of inh. population spiking to extracellular potassium accumulation
    modelp.gamma_a = 1.5;               % gamma_a: weight of astrocyte potassium clearance from the extracellular matrix
         
    % Parameters of KCl stimulus input    
    modelp.eta1 = 3;    % amplitude eta1 of the KCl stimulus input I
    modelp.eta2 = 0.4;  % sharpness eta2 of the nonlinearity (in space) of KCl stimulus input I
    modelp.tKCl = 2;    % initial time tKCl of KCl stimulus input I

    % Diffusion constant
    modelp.sigma = 1;

    %% Simulation parameters

    simp.nx = 2^11;                      % number of discretization nodes of the spatial grid
    simp.Lx = 100;                       % length of the spatial grid 
    simp.waveFrontThreshold = 0.91;      % threshold for CSD propagation detection in ve: we use the position where ve=wavefrontThreshold to find the propagation speed
    simp.errorVal = 0.2;                 % tolerance value if there is no point where ve=wavefrontThreshold. In this case, we pick up the point which provides the closest value to wavefrontThreshold within the tolerance determined by errorVal. 
    simp.T = 50;                         % final time of the simulation.
    simp.tStep = 1;                      % Output is sampled with tStep distances in time. Output is extracted at the time samples separated by tStep from each other. This does not determine the time integration step size, which is optimized automatically by ode45.
    simp.Cx = 13.6;                      % Scale factor of the speed. It is employed to fit the model to the experimental data in terms of propagation speed.
  
%% Save the parameters to .mat files modelp_control.mat and simp_control.mat
    
% Save the structure to a .mat file
save('Parameters/modelParameters/modelp_control.mat', '-struct', 'modelp');
save('Parameters/simulationParameters/simp_control.mat', '-struct', 'simp');

