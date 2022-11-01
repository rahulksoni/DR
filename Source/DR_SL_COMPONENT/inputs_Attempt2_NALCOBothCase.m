%% FileName = inputs_Attempt2_NALCOBothCase.m

%% Inputs for SQP program for Data Reconciliation
xlfilename          = '../../../NALCOProjectExcelDoNotRelocate';
SLSheetName    = 'SLInput';
M_F_0_range         = 'c2:c28';
U_F_0_range         = 'c29:c41';
M_SD_range          = 'b2:b28';
x_writerange        = 'd2:d41';
equality_write_range= 'h2:h22';

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
V_Zn_case       = 0;                            % 0 for one element case 
V_case          = 0;                            % 0 for one element case
Zn_case         = 0;                            % 0 for one element case

%% Convergence criteria
zero_tol          = 1e-06;                      % used to identify zero rows in Ru or inv(R1)*R2, and to detect duplicate initial values of variables
convergence_tol_M = 1e-05;                      % Convergence tolerance for measured variables
convergence_tol_U = 1e-05;                      % Convergence tolerance for unmeasured variables
convergence_flag  = 0;                          % flag turns to 1, if things get converged

%% Declaration of all possible symbolic variables based on total number of streams
n_stream = 38;                                  % total no. of streams in the flowsheet
for k=1:n_stream
syms(sprintf('Qs%d', k))                        % Solid flow rate
syms(sprintf('Ql%d', k))                        % Liquid flow rate
syms(sprintf('Cw%d', k))                        % Solid weight fraction in flow
syms(sprintf('CsV%d', k))                       % Solid V content
syms(sprintf('ClV%d', k))                       % Liquid V content
syms(sprintf('CsZn%d', k))                       % Solid Zn content
syms(sprintf('ClZn%d', k))                       % Liquid Zn content
end

%% Variables, modification and scaling (only if required)
scaling_factor = 1;

%% Variable values, standard deviations
% solid flow rate: measured and unmeasured variables 
M_Qs    =   [Qs1        Qs14        Qs16        Qs30        Qs33        Qs37];
% M_Qs_0  =  	[198.75		6.813752542	12.27228358	0.066666667	447.9618]*scaling_factor;                  
% SD_Qs   =   [6.57928849	2.043955902	2.93148299	0.026666667	0.30407903]*scaling_factor;       
% 
U_Qs    =   [Qs7        Qs12        Qs13        Qs15        Qs18        Qs26        Qs31        Qs36];
% U_Qs_0  =  	[198.8      110.8       87.8        103         120         150.4       150.3       150.1       150.2]*scaling_factor; 
% 
% % liquid flow rate: measured and unmeasured variables 
M_Ql    =   [Ql3        Ql4         Ql8         Ql9         Ql11        Ql17        Ql20        Ql28        Ql29        Ql32        Ql35];
% M_Ql_0  =  	[13.25      638.495		162.19875	62.1775		45.9		0.661666667	43.5		40          3868.185625	462.5		31.66666667]*scaling_factor;                  
% SD_Ql   =   [1.258305739 21.95329208 10.4770104	12.23141141	0.0001		0.099051053	1.732050808	28.28427125	552.0002724	94.64847243	5.773502692]*scaling_factor;       
% 
U_Ql    =   [Ql19      	Ql21        Ql27        Ql34        Ql38];
% U_Ql_0  =  	[10.0       50.0        50.1        4.0         2000.0]*scaling_factor;  
% 
% % solid fraction in flow: measured and unmeasured variables 
M_Cw    =   [Cw7        Cw12        Cw13        Cw14        Cw15        Cw16        Cw18        Cw26        Cw31        Cw33];
% M_Cw_0  =   [0.2226     0.2188      0.2772      0.081802    0.209       0.081801    0.2065      0.4234      0.2831      0.5902];
% SD_Cw   =   0.1*M_Cw_0;
% 
U_Cw    =   [];
% U_Cw_0  =   [];  
% % V content in solid
M_CsV    =   [CsV1       CsV7       CsV12       CsV13       CsV14       CsV15       CsV16       CsV18       CsV26       CsV31       CsV33       CsV36       CsV37];
M_CsV_0  =   [0.00521    0.0081     0.0097      0.0066      0.00092     0.0194      0.00091     0.00523     0.0013      0.0062      0.0074      0.0017      0.0008];     	
SD_CsV   =   0.1*M_CsV_0;     	

U_CsV    =   [];
U_CsV_0  =   []; 

% V content in liquid
M_ClV    =   [ClV4       ClV7        ClV8       ClV9        ClV11       ClV12       ClV13       ClV14       ClV15       ClV16       ClV18       ClV19       ClV20       ClV21       ClV26       ClV27       ClV29       ClV31       ClV33       ClV34];
M_ClV_0  =   [0.1247     0.0492      0.12675    0.12471     0.0297      0.0466      0.0303      0.1873      0.0636      0.18731     0.0556      0.0933      0.011       0.1077      0.0644      0.1261      0.0804      0.0371      0.0279      0.1257];  
SD_ClV   =   0.1*M_ClV_0;  

U_ClV    =   [];
U_ClV_0  =   [];

% Zn content in solid
M_CsZn   =   [CsZn1      CsZn7       CsZn12     CsZn13      CsZn14      CsZn15      CsZn16      CsZn18      CsZn26      CsZn31      CsZn33      CsZn36      CsZn37];
M_CsZn_0 =   [0.00521    0.0081      0.0097     0.0066      0.00092     0.0194      0.00091     0.00523     0.0013      0.0062      0.0074      0.0017      0.0008];     	
SD_CsZn  =   0.1*M_CsZn_0;     	

U_CsZn   =   [];
U_CsZn_0 =   []; 

% Zn content in liquid
M_ClZn   =   [ClZn4      ClZn7       ClZn8       ClZn9       ClZn11      ClZn12      ClZn13      ClZn14      ClZn15      ClZn16      ClZn18      ClZn19      ClZn20      ClZn21      ClZn26      ClZn27      ClZn29      ClZn31      ClZn33      ClZn34];
M_ClZn_0 =   [0.1247     0.0492      0.12675     0.12471     0.0297      0.0466      0.0303      0.1873      0.0636      0.18731     0.0556      0.0933      0.011       0.1077      0.0644      0.1261      0.0804      0.0371      0.0279      0.1257];  
SD_ClZn  =   0.1*M_ClZn_0;  

U_ClZn   =   [];
U_ClZn_0 =   [];

%% Extracting measured, unmeasured variables and standard deviations from excel file
M_F_0    = xlsread(xlfilename,SLSheetName,M_F_0_range)';   % extracting initial values of measured variables from excel file
U_F_0    = xlsread(xlfilename,SLSheetName,U_F_0_range)';   % extracting initial values of unmeasured variables from excel file
M_SD_F   = xlsread(xlfilename,SLSheetName,M_SD_range)';    % extracting initial values of measured variables from excel file

% Combining all M_sym and U_sym
if (FlowOnly_case == 1)
            M_sym   =  	[M_Qs       M_Ql        M_Cw];
            M_0     =	M_F_0;
            SD      =   M_SD_F;

            U_sym   =  	[U_Qs       U_Ql        U_Cw];
            U_0     =	U_F_0;
elseif (V_Zn_case == 1)
            M_sym   =  	[M_Qs       M_Ql        M_Cw       	M_CsV         M_ClV        M_CsZn        M_ClZn];
            M_0     =	[M_Qs_0     M_Ql_0      M_Cw_0      M_CsV_0       M_ClV_0      M_CsZn_0      M_ClZn_0];
            SD      =   [SD_Qs      SD_Ql       SD_Cw       SD_CsV        SD_ClV       SD_CsZn       SD_ClZn];

            U_sym   =  	[U_Qs       U_Ql        U_Cw       	U_CsV         U_ClV        U_CsZn        U_ClZn];
            U_0     =	[U_Qs_0     U_Ql_0      U_Cw_0      U_CsV_0       U_ClV_0      U_CsZn_0      U_ClZn_0];
elseif (V_case == 1)
            M_sym   =  	[M_Qs       M_Ql        M_Cw       	M_CsV         M_ClV];
            M_0     =	[M_Qs_0     M_Ql_0      M_Cw_0      M_CsV_0       M_ClV_0];
            SD      =   [SD_Qs      SD_Ql       SD_Cw       SD_CsV        SD_ClV];

            U_sym   =  	[U_Qs       U_Ql        U_Cw       	U_CsV         U_ClV];
            U_0     =	[U_Qs_0     U_Ql_0      U_Cw_0      U_CsV_0       U_ClV_0];
elseif (Zn_case == 1)
            M_sym   =  	[M_Qs       M_Ql        M_Cw       	M_CsZn        M_ClZn];
            M_0     =	[M_Qs_0     M_Ql_0      M_Cw_0      M_CsZn_0      M_ClZn_0];
            SD      =   [SD_Qs      SD_Ql       SD_Cw       SD_CsZn       SD_ClZn];

            U_sym   =  	[U_Qs       U_Ql        U_Cw       	U_CsZn        U_ClZn];
            U_0     =	[U_Qs_0     U_Ql_0      U_Cw_0      U_CsZn_0      U_ClZn_0];
else
            disp('Error: The case type selection such as Vanadium, Zinc or both had some issue')
            return
end

%% Provided answer
if (test_case_flag == 1)
    given_measured_variables_answer = given;
    length_given_measured_variables_answer = size(given_measured_variables_answer,2);
end

%% Commmon equality constraints
Common_equalities = [Qs1-Qs7;  % solid flow rate balancing
    Qs12+Qs14-Qs15;
    Qs15+Qs16-Qs18;
    Qs14+Qs16-Qs18-Qs30+Qs31;
    Qs13+Qs31-Qs33;
    Qs26-Qs36;
    Qs36-Qs37;
    Qs1+Qs30-Qs26-Qs33;

    % Solution mass balance (Balancing the liquid only flow rates)
    Ql3+Ql4-Qs7*(1/Cw7-1)+Ql8;
    Qs12*(1/Cw12-1)+Qs14*(1/Cw14-1)-Qs15*(1/Cw15-1);
    Qs15*(1/Cw15-1)+Qs16*(1/Cw16-1)+Ql17-Qs18*(1/Cw18-1)-Ql19;
    Qs14*(1/Cw14-1)+Qs16*(1/Cw16-1)-Qs18*(1/Cw18-1)+Qs31*(1/Cw31-1)-Ql34;
    Ql11-Qs13*(1/Cw13-1)-Qs31*(1/Cw31-1)-Ql32+Qs33*(1/Cw33-1);
    Ql19+Ql20-Ql21;
    Ql4+Ql9-Ql29+Ql38;
    (Ql27-Ql8)+Ql28-Ql29;
    Qs26*(1/Cw26-1)-Ql34+Ql35;
 
    % Total balance for digestion-sand separation and precipitation circuit
    Qs7/Cw7+Ql9+Ql11-Qs12/Cw12-Qs13/Cw13;
    Ql21-Qs26-Qs26*(1/Cw26-1)-Ql27;

    % Overall solid-liquid balance
    Ql3+Ql28+Ql17+Ql20+Ql32+Ql35-Qs33*(1/Cw33-1)-Ql38;

    % Additional supporting constraints
    Cw14-Cw16];

%% Valuable Vanadium content balance
V_equalities = [Qs1*CsV1+Ql4*ClV4-Qs7*CsV7-Qs7*(1/Cw7-1)*ClV7+Ql8*ClV8;    
    Qs7*CsV7+Qs7*(1/Cw7-1)*ClV7+Ql9*ClV9+Ql11*ClV11-Qs12*CsV12-Qs12*(1/Cw12-1)*ClV12-Qs13*CsV13-Qs13*(1/Cw13-1)*ClV13;    
    Qs12*CsV12+Qs12*(1/Cw12-1)*ClV12+Qs14*CsV14+Qs14*(1/Cw14-1)*ClV14-Qs15*CsV15-Qs15*(1/Cw15-1)*ClV15;
    Qs15*CsV15+Qs15*(1/Cw15-1)*ClV15+Qs16*CsV16+Qs16*(1/Cw16-1)*ClV16-Qs18*CsV18-Qs18*(1/Cw18-1)*ClV18-Ql19*ClV19;
    Ql19*ClV19+Ql20*ClV20-Ql21*ClV21;
    (Ql27-Ql8)*ClV8-Ql29*ClV29;
    Ql4*ClV4+Ql9*ClV9-Ql29*ClV29;
    Qs14*CsV14+Qs14*(1/Cw14-1)*ClV14+Qs16*CsV16+Qs16*(1/Cw16-1)*ClV16-Qs18*CsV18-Qs18*(1/Cw18-1)*ClV18+Qs31*CsV31+Qs31*(1/Cw31-1)*ClV31-Ql34*ClV34;
    Ql11*ClV11-Qs13*CsV13-Qs13*(1/Cw13-1)*ClV13-Qs31*CsV31-Qs31*(1/Cw31-1)*ClV31+Qs33*CsV33+Qs33*(1/Cw33-1)*ClV33;
    Qs26*CsV26+Qs26*(1/Cw26-1)*ClV26-Ql34*ClV34-Qs36*CsV36;
    Qs36*CsV36-Qs37*CsV37;
    Ql21*ClV21-Qs26*CsV26-Qs26*(1/Cw26-1)*ClV26-Ql27*ClV27;
    
    % Additional supporting constraints
    ClV4-ClV9;
    ClV8-ClV27;
    CsV14-CsV16;
    ClV14-ClV16; 
    CsV36-CsV37];
                    
%% Valuable Zinc content balance
Zn_equalities = [Qs1*CsZn1+Ql4*ClZn4-Qs7*CsZn7-Qs7*(1/Cw7-1)*ClZn7+Ql8*ClZn8;
    Qs7*CsZn7+Qs7*(1/Cw7-1)*ClZn7+Ql9*ClZn9+Ql11*ClZn11-Qs12*CsZn12-Qs12*(1/Cw12-1)*ClZn12-Qs13*CsZn13-Qs13*(1/Cw13-1)*ClZn13;    
    Qs12*CsZn12+Qs12*(1/Cw12-1)*ClZn12+Qs14*CsZn14+Qs14*(1/Cw14-1)*ClZn14-Qs15*CsZn15-Qs15*(1/Cw15-1)*ClZn15;
    Qs15*CsZn15+Qs15*(1/Cw15-1)*ClZn15+Qs16*CsZn16+Qs16*(1/Cw16-1)*ClZn16-Qs18*CsZn18-Qs18*(1/Cw18-1)*ClZn18-Ql19*ClZn19;
    Ql19*ClZn19+Ql20*ClZn20-Ql21*ClZn21;
    (Ql27-Ql8)*ClZn8-Ql29*ClZn29;
    Ql4*ClZn4+Ql9*ClZn9-Ql29*ClZn29;
    Qs14*CsZn14+Qs14*(1/Cw14-1)*ClZn14+Qs16*CsZn16+Qs16*(1/Cw16-1)*ClZn16-Qs18*CsZn18-Qs18*(1/Cw18-1)*ClZn18+Qs31*CsZn31+Qs31*(1/Cw31-1)*ClZn31-Ql34*ClZn34;
    Ql11*ClZn11-Qs13*CsZn13-Qs13*(1/Cw13-1)*ClZn13-Qs31*CsZn31-Qs31*(1/Cw31-1)*ClZn31+Qs33*CsZn33+Qs33*(1/Cw33-1)*ClZn33;
    Qs26*CsZn26+Qs26*(1/Cw26-1)*ClZn26-Ql34*ClZn34-Qs36*CsZn36;
    Qs36*CsZn36-Qs37*CsZn37;
    Ql21*ClZn21-Qs26*CsZn26-Qs26*(1/Cw26-1)*ClZn26-Ql27*ClZn27;
    
    % Additional supporting constraints
    ClZn4-ClZn9;
    ClZn8-ClZn27;
    CsZn14-CsZn16;
    ClZn14-ClZn16;
    CsZn36-CsZn37];

%% Equality constraints
if (FlowOnly_case == 1)
    equality_constraints = [Common_equalities];
elseif (V_Zn_case == 1)
    equality_constraints = [Common_equalities;
                            V_equalities;
                            Zn_equalities];
elseif (V_case == 1)
    equality_constraints = [Common_equalities;
                            V_equalities];
elseif (Zn_case == 1)
    equality_constraints = [Common_equalities;
                            Zn_equalities];
else
    disp('Error: The equalities due to case type selection such as Vanadium, Zinc or both had some issue')
    return
end

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