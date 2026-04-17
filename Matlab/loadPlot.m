
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This m.file is used to run the simulations of the astro-field model by
% using the "simulate" function, where we change one parameter in the
% course of a single time integration. The result corresponding to each
% parametere value is plotted as one block in the space-time plots.
%
% NOTE: THE NAME OF THE .MAT FILE TO BE LOADED SHOULD BE ENTERED IN
% MANUALLY BY THE USER. THE .MAT FILES ARE IN "DATA" FOLDER!!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;


% Add dependencies
addpath('Parameters/modelParameters');
addpath('Parameters/simulationParameters');
addpath('PlotFunctions');
addpath('SimulationFunctions');

%% Plot the saved simulation dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Space-time plots: 
% load('Data/fullSweep_gamma_a_from_10.5_to_10.5_changes_every_t_40_with_parStep_NaN.mat');  % UPDATE THE .MAT FILE NAME!!! loads x, tAll, all, idx
% 
% figure;
% PlotHistory(x, tAll, all, gcf, idx(:,1:4));
% 
% 
% figure;
% PlotHistoryTF(x, tAll, all, gcf, idx(:,5:7));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Propagation speed plots: 

% Load the data
% data = readmatrix('Data/GBZ_or_ISO_20260416_104348.csv');          % For ISO or GBZ case   
data = readmatrix('Data/ChR2_20260416_085443.csv');          % For ChR2 case   
% data = readmatrix('Data/GBZ_Plus_ISO_20260415_171644.csv');  % For GBZ + ISO case


% Extract columns
x = data(:,1);
y = data(:,2);

% Create figure
figure;

% Plot: blue markers ONLY (no line)
plot(x, y, 'b.', 'MarkerSize', 40);

% Labels (match your style)
xlabel('\gamma_{el} / \gamma_{il}');
ylabel('Speed');

% Axis limits (optional, based on your plot)
xlim([4 16]);
ylim([0 6]);

% Turn OFF grid explicitly
grid off;

% Optional: cleaner look like your figure
box on;
set(gca, 'FontSize', 12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
