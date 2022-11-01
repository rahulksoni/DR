%% FileName = inputs_Attempt2_NALCOBothCase.m

%% Inputs for SQP program for Data Reconciliation
% xlfilename          = '../../../NALCOProjectExcelDoNotRelocate';
% SLSheetName         = 'SLInput';
% M_F_0_range         = 'c2:c28';
% U_F_0_range         = 'c29:c41';
% M_SD_range          = 'b2:b28';
% x_writerange        = 'd2:d41';
% equality_write_range= 'h2:h22';

%% Declaration and value of all the inputs need to be done here
% File in which the content of this file need to be copied is:
% inputs.m

%% Instructions and note points
% Capital M and N refers to measured and unmeasured variables, the small
% notications as m and n may refer to matrix sizes

%% Case Dependent Comment: 
% 1. The present case is based on paper titled Nonlinear data reconciliation in gold processing plants authored by L.R.P. de Andrade Lima.
% 2. The total number of streams are 38
% 3. There are following types of variable in the case
% 	Solid flow rate: Qs
% 	Liquid flow rate: Ql
% 	Solid fraction in slurry: Cw
% 	Solid element content: Cs
% 	Liquid element content: Cl
% 4. So, there are 5 types of variable for which measured and unmeasured values are available 

%% Specific to file comment


%% If this is a test case with provided answer 
test_case_flag = 0;                             % otherwise, put 1

%% If the case is for only Vanadium or only Zinc or both Vanadium and Zinc (make either of these three one)
FlowOnly_case   = 1;                            % to only reconclide flows without element content balance

%% Convergence criteria
zero_tol          = 1e-06;                      % used to identify zero rows in Ru or inv(R1)*R2, and to detect duplicate initial values of variables
convergence_tol_M = 1e-05;                      % Convergence tolerance for measured variables
convergence_tol_U = 1e-05;                      % Convergence tolerance for unmeasured variables
convergence_flag  = 0;                          % flag turns to 1, if things get converged

%% Declaration of all possible symbolic variables based on total number of streams
n_stream = 2;                                  % total no. of streams in the flowsheet
for k=1:n_stream
syms(sprintf('x%d', k))                        % Solid flow rate
% syms(sprintf('Ql%d', k))                        % Liquid flow rate
% syms(sprintf('Cw%d', k))                        % Solid weight fraction in flow
% syms(sprintf('CsV%d', k))                       % Solid V content
% syms(sprintf('ClV%d', k))                       % Liquid V content
% syms(sprintf('CsZn%d', k))                       % Solid Zn content
% syms(sprintf('ClZn%d', k))                       % Liquid Zn content
end

%% Variables, modification and scaling (only if required)
scaling_factor = 1;

%% Variable values, standard deviations
% solid flow rate: measured and unmeasured variables 
M_x    =   [x1 x2];
M_x_0  =   [3 6]*scaling_factor;                  
SD_x   =   M_x_0*0.5;       

% Combining all M_sym and U_sym

            M_sym   =  	M_x;
            U_sym   =  [];
            U_0     = [];
            M_0     =	M_x_0;
            SD      =   SD_x;

%% Provided answer
if (test_case_flag == 1)
    given_measured_variables_answer = given;
    length_given_measured_variables_answer = size(given_measured_variables_answer,2);
end

%% Commmon equality constraints
   equality_constraints = [x1*x1-x2];


%% All inequality constraint need to be in =< form. e.g. a+b =< 0
% Note that successive linearization cannot handle inequality constraint
inequality_constraints = [-M_sym';      % Positive value constraint
        -U_sym';

        % Bound constraints
        []'-100;
    
        % Inequality physical constraints
        ];
    
%% Mass balancing equations for calculation of non-redundant variables
non_redundants_equations = [%U_sym(1)+M_sym(2);             % mass balancing equation in form of only equating right part, e.g. equation F1 = F2+F3 shall be written as only F2+F3
                            %M_sym(2)
                            ];  