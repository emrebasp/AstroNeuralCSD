
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% RUN SIMULATIONS -- GBZ + ChR2 case
%
% This m.file is used to run the simulations of the astro-field model by
% using the "simulate" function, where we change one parameter in the
% course of a single time integration. The result corresponding to each
% parametere value is plotted as one block in the space-time plots.
%
% This code was designed and edited by
% Emre Baspinar, Inria, MathNeuro Team, France,
% emre.baspinar@inria.fr,
% Daniele Avitabile, VU Amsterdam, Netherlands
% d.avitabile@vu.nl .
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   GLOSSARY FOR THE PARAMETERS WITH THEIR CONTROL VALUES
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Parameters of population transfer functions Se, Si, Sa    
% modelp.theta1_e = 1.7;              % theta1_e: initial potassium level of hyperactivation of the exc. population
% modelp.theta2_e = 1.8;              % theta2_e: initial potassium level of depolarization block of the exc. population
% modelp.theta1_i = 1.2;              % theta1_i: initial potassium level of hyperactivation of the inh. population
% modelp.theta2_i = 1.4;              % theta2_i: initial potassium level of depolarization block of the inh. population
% modelp.theta1_a = 1.0;              % theta1_a: initial potassium level of increased potassium clearance of the astro. population
% modelp.theta2_a = 1.2;              % theta2_a: initial potassium level of no potassium clearance  
% modelp.theta3_e  = 100;             % theta3_e: sharpness of the nonlinearity of S_e 
% modelp.theta3_i  = 10;              % theta3_i: sharpness of the nonlinearity of S_i 
% modelp.theta3_a  = 100;             % theta3_a: sharpness of the nonlinearity of S_a 
% modelp.theta4  = 0.3;               % theta4: threshold of the transfer functions S_e, S_i, S_a
% 
% % Connection weights
% modelp.gamma_ei = 4;                % gamma_ei: connection weight from inh. to exc. neurons
% modelp.gamma_ie = 0.2;              % gamma_ie: connection weight from exc. to inh. neurons
% modelp.gamma_ee = 1.0;              % gamma_ee: connection weight from exc. to exc. neurons
% modelp.gamma_ii = 4;                % gamma_ii: connection weight from inh. to inh. neurons
% modelp.gamma_an = 0.5;              % gamma_an: connection weight from all neurons to astrocytes
% modelp.gamma_aa = 0.5;              % gamma_aa: connection weight from astrocytes to astrocytes
% 
% % Parameters of potassium-dependent activation functions Fe, Fi, Fa
% modelp.alpha1_e  = 100;             % alpha1_e: sharpness of the nonlinerarity of Fe in ve equation 
% modelp.alpha1_i  = 10;              % alpha1_i: sharpness of the nonlinerarity of Fi in vi equation
% modelp.alpha1_a  = 100;             % alpha1_a: sharpness of the nonlinerarity of Fa in ra equation
% modelp.alpha2_e = 1.05;             % alpha2_e: potassium threshold of Fe function appearing in ve equation
% modelp.alpha2_i = 1.05;             % alpha2_i: potassium threshold of Fi function appearing in vi equation
% modelp.alpha2_a = 1.05;             % alpha2_a: potassium threshold of Fa function appearing in ra equation
% 
% % Weights of potassium-dependent activation functions Fe, Fi, Fa
% modelp.gamma_eK = 1;         % weight for Fe
% modelp.gamma_iK = 1;         % weight for Fi
% modelp.gamma_aK = 1;         % weight for Fa
% 
% % Parameters of potassium source function rho
% modelp.beta1_e = 10;                % beta1_e: sharpness of the nonlinearity of the rho function weighted by gamma_e in potassium equation
% modelp.beta2_e = 1.5;               % beta2_e: threshold value for the rho function weighted by gamma_e in potassium equation
% modelp.beta1_i = 10;                % beta1_i: sharpness of the nonlinearity of the rho function weighted by gamma_i in potassium equation
% modelp.beta2_i = 1.5*0.9975;        % beta2_i: threshold value for the rho function weighted by gamma_i in potassium equation
% modelp.beta1_a = 10;                % beta1_i: sharpness of the nonlinearity of the rho function weighted by gamma_a in potassium equation
% modelp.beta2_a = 1.5;               % beta2_a: threshold value for the rho function weighted by c3 in potassium equation
% 
% % Weights of potassium source function rho
% modelp.gamma_e = 18.6;              % gamma_e: weight of contribution of exc. population spiking to extracellular potassium accumulation
% modelp.gamma_i = 18.6;              % gamma_i: weight of contribution of inh. population spiking to extracellular potassium accumulation
% modelp.gamma_a = 1.5;               % gamma_a: weight of astrocyte potassium clearance from the extracellular matrix
% 
% % Parameters of KCl stimulus input    
% modelp.eta1 = 3;    % amplitude eta1 of the KCl stimulus input I
% modelp.eta2 = 0.4;  % sharpness eta2 of the nonlinearity (in space) of KCl stimulus input I
% modelp.tKCl = 2;    % initial time tKCl of KCl stimulus input I
% 
% % Diffusion constant
% modelp.sigma = 1;
% 
% %% Simulation parameters
% 
% simp.nx = 2^11;                      % number of discretization nodes of the spatial grid
% simp.Lx = 100;                       % length of the spatial grid 
% simp.waveFrontThreshold = 0.91;      % threshold for CSD propagation detection in ve: we use the position where ve=wavefrontThreshold to find the propagation speed
% simp.errorVal = 0.2;                 % tolerance value if there is no point where ve=wavefrontThreshold. In this case, we pick up the point which provides the closest value to wavefrontThreshold within the tolerance determined by errorVal. 
% simp.T = 50;                         % final time of the simulation.
% simp.tStep = 1;                      % Output is sampled with tStep distances in time. Output is extracted at the time samples separated by tStep from each other. This does not determine the time integration step size, which is optimized automatically by ode45.
% simp.Cx = 13.6;                      % Scale factor of the speed. It is employed to fit the model to the experimental data in terms of propagation speed.

%==========================================================================


clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Initialization
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except parIndex

warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle') % Turn off this warning, it does not have any influence on the code functionality.

% Add dependencies
addpath('Parameters/modelParameters');
addpath('Parameters/simulationParameters');
addpath('SimulationFunctions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the model parameters
modelp = load('modelp_control.mat'); 

%==========================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the simulation parameters
simp = load('simp_control.mat');
simp.Lx = simp.Lx;
simp.nx = simp.nx;

%==========================================================================


%% Space-time plots

% Grid parameters
Lx = simp.Lx;
nx = simp.nx;
hx = Lx/(nx-1);             % distance between two discretization nodes 
x = -Lx/2 + [0:nx-1]'*hx;   % discretization nodes
T = simp.T;
tStep = simp.tStep;
tSpan = [0:tStep:T];        % Output is sampled at the instants in tSpan.

% Indexing the state variables
iVe = [1:nx]'; % Indexes for ve vector (exc. population)
iVi = nx+iVe;  % Indexes for vi vector (inh. population)
iRa = nx+iVi;  % Indexes for ra vector (astro. population)
iCK = nx+iRa;  % Indexes for cK vector (extracellular potassium concentration)

% All indexes in idx, which will have size(idx)=(nx, 7).
idx = zeros(nx,7);

idx(:,1) = iVe;                % ve
idx(:,2) = iVi;                % vi
idx(:,3) = iRa;                % ra
idx(:,4) = iCK;                % cK

idx(:,5) = 4*nx + (1:nx)';     % Se block in combined data vector
idx(:,6) = 5*nx + (1:nx)';     % Si block
idx(:,7) = 6*nx + (1:nx)';     % Sa block

% Parameter ranges for varied connection weights gamma_ei, gamma_ii and
% varied gamma_i
% nOfSamples = 61;
% initRange1 = 2; % for gamma_ei, gamma_ii
% finRange1 = 12;
% initRange2 = 4; % for gamma_i
% finRange2 = 34;
nOfSamples = 61;
initRange1 = 1; % for gamma_ei, gamma_ii
finRange1 = 11;
initRange2 = 12; % for gamma_i
finRange2 = 102;

parRange1 = linspace(initRange1,finRange1,nOfSamples); % for gamma_ei, gamma_ii
parRange1 = flip(parRange1); % reverse the range since gamma_ei, gamma_ii will be decreasing 
parRange2 = linspace(initRange2,finRange2,nOfSamples); % for gamma_i

% Initial conditions
ve0 = zeros(size(x));
vi0 = zeros(size(x));
ra0 = zeros(size(x));
cK0  = zeros(size(x));
z0  = [ve0; vi0; ra0; cK0];

modeSpeed = 1; % Compute the propagation speed in "simulate" function.
propagationSpeed=zeros(1, nOfSamples); % Initialize a vector to save the CSD propagation speed values

% Perform the simulation for each gamma_ei, gamma_cii and for each gamma_i
for i=1:nOfSamples
    
    % Keep track of the simulations
    fprintf('Iteration %d out of %d\n', i, nOfSamples);

    % Vary gamma_ei, gamma_ii and gamma_i
    modelp.gamma_ei = parRange1(i);
    modelp.gamma_ii = parRange1(i);
    modelp.gamma_i  = parRange2(i);

    % % Keep track of the value of gamma_ei,gamma_ii and gamma_i
    % fprintf('gamma_{ei}=gamma_{ii]=%.3f and gamma_i=%.3f\n', modelp.gamma_ei, modelp.gamma_i);


    result = simulate(modelp, simp, tSpan, z0, modeSpeed);     
    propagationSpeed(i) = result.propagationSpeed; % save the speed

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Plot and save the propagation speeds
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the propagation speed wrt to varied gamma_ei, gamma_ii values
figure;
plot(parRange1, propagationSpeed, 'bo','MarkerFaceColor', 'b', 'LineWidth', 2);
axis([parRange1(end) parRange1(1) 0 6])
set(gca, 'XDir', 'reverse'); % Reverse the X-axis since gamma_ei, gamma_ii are decreasing
grid off;
xlabel('\gamma_i');
ylabel('Speed');

% Uncomment below to save the figure
% Define the folder name where the figure will be saved
folderName = 'Plots'

% Get the directory of the currently running .m file
scriptDir = fileparts(mfilename('fullpath'));

% Create the folder if it doesn't exist
savePath = fullfile(scriptDir, folderName);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

% Generate the timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS'); % Format: YYYYMMDD_HHMMSS

% Define the filename with the timestamp
fileName = sprintf('GBZ_Plus_ISO_%s.png', timestamp);

% Full path to save the file
fullFilePath = fullfile(savePath, fileName);

% Save the figure as PNG
saveas(gcf, fullFilePath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the folder name where the CSV file will be saved
folderName = 'Data';

% Get the directory of the currently running .m file
scriptDir = fileparts(mfilename('fullpath'));

% Create the folder if it doesn't exist
savePath = fullfile(scriptDir, folderName);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

% Generate the timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS'); % Format: YYYYMMDD_HHMMSS

% Define the filename with the timestamp
fileName = sprintf('GBZ_Plus_ISO_%s.csv', timestamp);

% Full path to save the file
fullFilePath = fullfile(savePath, fileName);

% Write data to the CSV file
data = [parRange1' propagationSpeed']; % Combine data into a matrix (column-wise)
writematrix(data, fullFilePath, 'WriteMode', 'append'); 

