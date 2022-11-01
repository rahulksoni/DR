function [c, ceq] = DR_SL_Constraints( x )
% Both inequality and equality contraints are called from this file
% Note: this file has no value for successive linearization scheme

%% Setting up global and local variables 
global M_sym U_sym length_M length_U result input_file_name
M_sym = M_sym; U_sym = U_sym; length_M = length_M; length_U = length_U; result = result; input_file_name = input_file_name;

%% Variable value assignment
n_rows_result = size(result,1);                             % the counting till last iteration
M = x(1:length_M);                                          % values of measured variables from last iteration 
U = result(n_rows_result,length_M+1:length_M+length_U);     % values of unmeasured variables from last iteration

%% Input data/material is copied here for the purpose of getting equality and inequality constraints
eval(input_file_name);

%% Filling up symbolic variables with original values for later subsitutions
for i = 1:length_M
    evalc(sprintf('%s=M(%d)',char(M_sym(i)),i));
end
for i = 1:length_U
    evalc(sprintf('%s=U(%d)',char(U_sym(i)),i));
end

%% equality and inequality constraints with their real values
c   = double(subs(inequality_constraints));
ceq = double(subs(equality_constraints));