function S = PopulationTransferFcn(v, cK, theta1, theta2, theta3, theta4)

% XI = POPULATIONTRANSFERFCN(V, CK, THETA1, THETA2, THETA3, THETA4) provides the
% output XI of the population transfer function of corresponding cell population.

% POPULATIONTRANSFERFCN(V, CK, THETA1, THETA2, THETA3, THETA4): V is the average
% membrane potential or the average gliotransmitter release activity. These are
% the inputs for the transfer functions of neural populations or for the
% astrocyte population. V can be one of the following inputs: ve, vi,
% ve+vi, ra. CK is the extracellular potassium concentration. THETA1 and THETA2 give the
% extracellular potassium levels which determine normal activity, high and
% very high activity phases of the transfer function. THETA3 determines 
% the sharpness of the nonlinearities. THETA4 is the input threshold. 
% It determines the input value corresponding to the inflection
% point of the nonlinearities of the transfer function.  
% Initialize a vector for the output values of the transfer function

S = zeros(size(v));

% Transfer function formula
S =  (cK<theta1) .* 0.5 ./ (1 + exp(- theta3 * (v - theta4) ))...
 +  (cK>=theta1) .* (cK<= theta2) ./ (1 + exp(- theta3 * (v - theta4) ));

end