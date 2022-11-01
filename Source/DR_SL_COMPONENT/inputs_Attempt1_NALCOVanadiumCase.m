%% FileName = inputs_Attempt1_NALCOVanadiumCase.m





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
syms(sprintf('Cs%d', k))                        % Solid Zn content
syms(sprintf('Cl%d', k))                        % Liquid Zn content
end

%% Variables, modification and scaling (only if required)

%% Si

%% Variable values, standard deviations
% solid fraction in flow: measured and unmeasured variables 
M_Cw    =   [Cw7        Cw10        Cw12        Cw13        Cw14        Cw15        Cw16        Cw18        Cw26        Cw31        Cw33];
M_Cw_0  =   [0.2226     0.3607      0.2188      0.2772      0.081802    0.209       0.081801    0.2065      0.4234      0.2831      0.5902];
SD_Cw   =   0.1*M_Cw_0;

U_Cw    =   [];
U_Cw_0  =   [];  

% solid flow rate: measured and unmeasured variables 
M_Qs    =   [Qs1        Qs14        Qs16        Qs30        Qs33];
M_Qs_0  =  	[149.583    4.668       9.449       0.0467      333.536];                  
SD_Qs   =   [99.929     3.531       6.917       0.04        236.325];       

U_Qs    =   [Qs7        Qs10        Qs12        Qs13        Qs15        Qs18        Qs26        Qs31        Qs36        Qs37];
U_Qs_0  =  	[149.58     149.581     70.0        80.0        135.0       150.5       150.4       150.3       150.1       150.2]; 

% liquid flow rate: measured and unmeasured variables 
M_Ql    =   [Ql3        Ql4         Ql8         Ql9         Ql11        Ql17        Ql20        Ql28        Ql29        Ql32        Ql35];
M_Ql_0  =  	[9.667      475.475     121.283     49.115      34.425      0.514       33.01       25.02       2860.442    362.5       24.167];                  
SD_Ql   =   [6.526      318.61      81.52       34.32       22.95       0.353       22.04       30          1982.27     256.17      17.07];       

U_Ql    =   [Ql19      	Ql21        Ql27        Ql34        Ql38];
U_Ql_0  =  	[10.0       50.0        50.1        4.0         2000.0];  

% V content in solid
M_Cs    =   [Cs1        Cs7         Cs10        Cs12        Cs13        Cs14        Cs15        Cs16        Cs18        Cs26        Cs31        Cs33        Cs36        Cs37];
M_Cs_0  =   [0.00521    0.0081      0.0281      0.0097      0.0066      0.00092     0.0194      0.00091     0.00523     0.0013      0.0062      0.0074      0.0017      0.0008];     	
SD_Cs   =   0.1*M_Cs_0;     	

U_Cs    =   [];
U_Cs_0  =   []; 

% V content in liquid
M_Cl    =   [Cl4      	Cl7         Cl8         Cl9         Cl10        Cl11        Cl12        Cl13        Cl14        Cl15        Cl16        Cl18        Cl19        Cl20        Cl21        Cl26        Cl27        Cl29        Cl31        Cl33        Cl34];
M_Cl_0  =   [0.1247     0.0492      0.12675     0.12471     0.0538      0.0297      0.0466      0.0303      0.1873      0.0636      0.18731     0.0556      0.0933      0.011       0.1077      0.0644      0.1261      0.0804      0.0371      0.0279      0.1257];  
SD_Cl   =   0.1*M_Cl_0;  

U_Cl   =    [];
U_Cl_0 =    [];

% Combining all M_sym and U_sym
M_sym   =  	[M_Qs       M_Ql        M_Cw       	M_Cs        M_Cl];
M_0     =	[M_Qs_0     M_Ql_0      M_Cw_0      M_Cs_0      M_Cl_0];
SD      =   [SD_Qs      SD_Ql       SD_Cw       SD_Cs       SD_Cl];

U_sym   =  	[U_Qs       U_Ql        U_Cw       	U_Cs        U_Cl];
U_0     =	[U_Qs_0     U_Ql_0      U_Cw_0      U_Cs_0      U_Cl_0];

%% Provided answer
if (test_case_flag == 1)
    given_measured_variables_answer = given;
    length_given_measured_variables_answer = size(given_measured_variables_answer,2);
end

%% Equality constraints function definition
equality_constraints = [Qs1-Qs7;                % solid flow rate balancing
                        Qs7-Qs10;
                        Qs10-Qs12-Qs13;
                        Qs12+Qs14-Qs15;
                        Qs15+Qs16-Qs18;
                        Qs14+Qs16-Qs18-Qs30+Qs31;
                        Qs13+Qs31-Qs33;
                        Qs26-Qs36;
                        Qs36-Qs37;
        
                        % Solution mass balance (Balancing the liquid only flow rates)
                        Ql3+Ql4-Qs7*(1/Cw7-1)+Ql8;
                        Qs7*(1/Cw7-1)+Ql9-Qs10*(1/Cw10-1);
                        Qs10*(1/Cw10-1)+Ql11-Qs12*(1/Cw12-1)-Qs13*(1/Cw13-1);
                        Qs12*(1/Cw12-1)+Qs14*(1/Cw14-1)-Qs15*(1/Cw15-1);
                        Qs15*(1/Cw15-1)+Qs16*(1/Cw16-1)+Ql17-Qs18*(1/Cw18-1)-Ql19;
                        Ql19+Ql20-Ql21;
                        (Ql27-Ql8)+Ql28-Ql29;%Ql8-Ql27-Ql28+Ql29;
                        Ql4+Ql9-Ql29+Ql38;
                        Qs14*(1/Cw14-1)+Qs16*(1/Cw16-1)-Qs18*(1/Cw18-1)+Qs31*(1/Cw31-1)-Ql34;
                        Ql11-Qs13*(1/Cw13-1)-Qs31*(1/Cw31-1)-Ql32+Qs33*(1/Cw33-1);
                        Qs26*(1/Cw26-1)-Ql34+Ql35;
                        
                        % Solid-liquid combined balance for precipitation
                        Ql21-Qs26-Qs26*(1/Cw26-1)-Ql27;
                        
                        % Valuable content balance
                        Qs1*Cs1+Ql4*Cl4-Qs7*Cs7-Qs7*(1/Cw7-1)*Cl7+Ql8*Cl8;
                        Qs7*Cs7+Qs7*(1/Cw7-1)*Cl7+Ql9*Cl9-Qs10*Cs10-Qs10*(1/Cw10-1)*Cl10;
                        Qs10*Cs10+Qs10*(1/Cw10-1)*Cl10+Ql11*Cl11-Qs12*Cs12-Qs12*(1/Cw12-1)*Cl12-Qs13*Cs13-Qs13*(1/Cw13-1)*Cl13;
                        Qs12*Cs12+Qs12*(1/Cw12-1)*Cl12+Qs14*Cs14+Qs14*(1/Cw14-1)*Cl14-Qs15*Cs15-Qs15*(1/Cw15-1)*Cl15;
                        Qs15*Cs15+Qs15*(1/Cw15-1)*Cl15+Qs16*Cs16+Qs16*(1/Cw16-1)*Cl16-Qs18*Cs18-Qs18*(1/Cw18-1)*Cl18-Ql19*Cl19;
                        Ql19*Cl19+Ql20*Cl20-Ql21*Cl21;
                        (Ql27-Ql8)*Cl8-Ql29*Cl29;
                        Ql4*Cl4+Ql9*Cl9-Ql29*Cl29;
                        Qs14*Cs14+Qs14*(1/Cw14-1)*Cl14+Qs16*Cs16+Qs16*(1/Cw16-1)*Cl16-Qs18*Cs18-Qs18*(1/Cw18-1)*Cl18+Qs31*Cs31+Qs31*(1/Cw31-1)*Cl31-Ql34*Cl34;
                        Ql11*Cl11-Qs13*Cs13-Qs13*(1/Cw13-1)*Cl13-Qs31*Cs31-Qs31*(1/Cw31-1)*Cl31+Qs33*Cs33+Qs33*(1/Cw33-1)*Cl33;
                        Qs26*Cs26+Qs26*(1/Cw26-1)*Cl26-Ql34*Cl34-Qs36*Cs36;
                        Qs36*Cs36-Qs37*Cs37;
                        Ql21*Cl21-Qs26*Cs26-Qs26*(1/Cw26-1)*Cl26-Ql27*Cl27;
                        
                        % Additional supporting constraints
                        Cl4-Cl9;
                        %Cl8-Cl29;
                        Cl8-Cl27;
                        Cs14-Cs16;
                        Cl14-Cl16;
                        Cw14-Cw16;
                        Cs26-Cs36;
                        Cs36-Cs37;
                        ];

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
                        
%% RESULTING COMMENTS ON SOLVING THE GOLD CASE