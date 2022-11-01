%% FileName = inputs_case_GoldCase.m





%% Declaration and value of all the inputs need to be done here
% File in which the content of this file need to be copied is:
% inputs.m

%% Instructions and note points
% Capital M and N refers to measured and unmeasured variables, the small
% notications as m and n may refer to matrix sizes

%% Case Dependent Comment: Following is the example based on case 1 in paper titled Nonlinear data reconciliation in gold processing plants by L.R.P. de Andrade Lima
% 1. The present case is based on paper titled Nonlinear data reconciliation in gold processing plants authored by L.R.P. de Andrade Lima.
% 2. The total number of streams are 23
% 3. There are following types of variable in the case
% 	Solid flow rate: Qs
% 	Liquid flow rate: Ql
% 	Solid concentration: Cw
% 	Solid gold content: Cs
% 	Liquid gold content: Cl
% 4. So, there are 5 types of variable for which measured and unmeasured values are available 
% 5. So, the variable vector size of 23x5=115

%% Specific to file
% major modifications are done here as the equations given in paper are not
% right, there are serious mismatch in the units. ALso, the standard
% deviation values were wrong and major culprit to the unsuccessful results
% from the case.

%% If this is a test case with provided answer 
test_case_flag = 1;                             % otherwise, put 0

%% Convergence criteria
zero_tol          = 1e-06;                      % used to identify zero rows in Ru or inv(R1)*R2, and to detect duplicate initial values of variables
convergence_tol_M = 1e-5;                       % Convergence tolerance for measured variables
convergence_tol_U = 1e-5;                       % Convergence tolerance for unmeasured variables
convergence_flag  = 0;                          % flag turns to 1, if things get converged

%% Declaration of all possible symbolic variables based on total number of streams
n_stream = 23;                                  % total no. of streams in the flowsheet
for k=1:n_stream
syms(sprintf('Qs%d', k))                        % Solid flow rate
syms(sprintf('Ql%d', k))                        % Liquid flow rate
syms(sprintf('Cw%d', k))                        % Solid fraction in flow
syms(sprintf('Cs%d', k))                        % Solid gold content
syms(sprintf('Cl%d', k))                        % Liquid gold content
syms(sprintf('Au%d', k))                        % Gold value 
end

%% Variables, modification and scaling
% scaling factors
Qs_scale_factor = 1000000/1000000;                                                  % converting t/h to gram/h; Modification 3 & 4 (scaling down the Qs, Ql and Au values by 1000000);
Cw_scale_factor = 1;                                                                % no change
Ql_scale_factor = 1000000/1000000;                                                  % converting t/h to gram/h; Modification 3 & 4 (scaling down the Qs, Ql and Au values by 1000000)
Cs_scale_factor = 1/1000000;                                                        % converting mg/kg to gram/gram 
Cl_scale_factor = 1/1000000;                                                        % converting mg/L to gram/gram (consider fluid density as of water because gold content is very low)
Au_scale_factor = 1/1000000;                                                        % no change; Modification 3 & 4 (scaling down the Qs, Ql and Au values by 1000000);

% standard deviation modification factors
Qs_SD_mod       = 0.2;
Cw_SD_mod       = 0.5;
Ql_SD_mod       = 0.35;
Cs_SD_mod       = 0.5;
Cl_SD_mod       = 0.5;
Au_SD_mod       = 0.5;

% solid flow rate
M_Qs    =   [Qs1];
M_Qs_0  =  	[55.981]*Qs_scale_factor;                  
SD_Qs   =   [11.196]*Qs_scale_factor*Qs_SD_mod;       
giv_Qs  =   [54.03]*Qs_scale_factor; 
% U_Qs    = [Qs2      	Qs4         Qs8         Qs10        Qs11        Qs13        Qs16];
% U_Qs_0  = [54.0301    54.0302     54.0303     54.0304     54.0305     54.0306     54.0307]*Qs_scale_factor; 
% Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
U_Qs    =   [Qs2      	Qs4                     Qs10                    Qs13        Qs16];
U_Qs_0  =  	[54.0301    54.0302                 54.0304                 54.0306     54.0307]*Qs_scale_factor; 

% solid fraction in flow
% M_Cw    = [Cw2        Cw4         Cw8         Cw10        Cw11        Cw13        Cw16]*Cw_scale_factor;
% M_Cw_0  = [0.29       0.5529      0.5265   	0.5237      0.5151      0.5231      0.5143]*Cw_scale_factor;
% SD_Cw   = [0.058      0.11058     0.1053      0.10474     0.10302     0.10462     0.10286]*Cw_scale_factor;
% giv_Cw  = [0.3709     0.5724      0.5724      0.5236      0.5236      0.5822      0.6218]*Cw_scale_factor;
% Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
M_Cw    =   [Cw2        Cw4                    	Cw10                    Cw13        Cw16];
M_Cw_0  =   [0.29       0.5529                 	0.5237              	0.5231      0.5143]*Cw_scale_factor;
SD_Cw   =   [0.058      0.11058                 0.10474                 0.10462     0.10286]*Cw_scale_factor*Cw_SD_mod;
giv_Cw  =   [0.3709     0.5724                  0.5236                  0.5822      0.6218]*Cw_scale_factor;
U_Cw    =   []*Cw_scale_factor;
U_Cw_0  =   []*Cw_scale_factor;

% liquid flow rate
% M_Ql    =   [Ql3      Ql5         Ql6         Ql7         Ql9         Ql12        Ql14        Ql19        Ql21];
% M_Ql_0  =   [139.5    28.59       154.94      121.19      29.17       145.91      144.22      24.72       154.92];
% SD_Ql   =   [97.65    20.013      108.458     84.833      20.419      102.137     100.954     17.304      108.444];
% giv_Ql  =   [91.63    41.49       167.52      74.78       50.27       115.53      105.14      40.75       167.52];
% Modification 1: Converting unmeasured variable Ql20 to measured variable
% M_Ql    =   [Ql3      Ql5         Ql6         Ql7         Ql9         Ql12        Ql14        Ql19        Ql21        Ql20];
% M_Ql_0  =   [139.5    28.59       154.94      121.19      29.17       145.91      144.22      24.72       154.92      20.9]*Ql_scale_factor;     
% SD_Ql   =   [97.65    20.013      108.458     84.833      20.419      102.137     100.954     17.304      108.444     0.05]*Ql_scale_factor;     
% giv_Ql  =   [91.63    41.49       167.52      74.78       50.27       115.53      105.14      40.75       167.52      20.97]*Ql_scale_factor;    
% Modification 2: Converting unmeasured variable Ql22 to measured variable
M_Ql    =   [Ql3        Ql5         Ql6         Ql7         Ql9         Ql12        Ql14        Ql19        Ql21        Ql20        Ql22];
M_Ql_0  =   [139.5     	28.59       154.94      121.19      29.17       145.91      144.22      24.72       154.92      20.9        0.0000001]*Ql_scale_factor;    
SD_Ql   =   [97.65     	20.013      108.458     84.833      20.419      102.137     100.954     17.304      108.444     0.05        0.0004]*Ql_scale_factor*Ql_SD_mod;  
giv_Ql  =   [91.63      41.49       167.52      74.78       50.27       115.53      105.14      40.75       167.52      20.97       0.0]*Ql_scale_factor;  
% U_Ql    =   [Ql15     Ql17       	Ql18        Ql20       	Ql22];
% U_Ql_0  =   [29.91    11.89       12.11       20.97       0.0];
% Modification 1: Converting unmeasured variable Ql20 to measured variable
% U_Ql    =   [Ql15     Ql17       	Ql18                   	Ql22];
% U_Ql_0  =   [29.91    11.89       12.11                   0.0]*Ql_scale_factor;   
% Modification 2: Converting unmeasured variable Ql22 to measured variable
U_Ql    =   [Ql15     	Ql17       	Ql18                   	];
U_Ql_0  =   [29.91      11.89       12.11                   ]*Ql_scale_factor;   

% gold content in solid
% M_Cs    =   [Cs1       	Cs2         Cs4         Cs10        Cs11        Cs13        Cs16];
% M_Cs_0  =   [16.459    	3.733       3.5036      1.293       0.5839      0.4171      0.3545]*Cs_scale_factor;     	
% SD_Cs   =   [14.8131   	2.6131      2.52252     0.9051      0.40873     0.29197     0.24815]*Cs_scale_factor;     	
% giv_Cs  =   [5.006      3.0724      2.7241      1.3986      0.6247      0.6247      0.1972]*Cs_scale_factor;    
% U_Cs    =   [Cs8];
% U_Cs_0  =   [1.3986]*Cs_scale_factor; 
% Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
M_Cs    =   [Cs1       	Cs2         Cs4         Cs10                    Cs13        Cs16];
M_Cs_0  =   [16.459    	3.733       3.5036      1.293                   0.4171      0.3545]*Cs_scale_factor;     	
SD_Cs   =   [14.8131   	2.6131      2.52252     0.9051                  0.29197     0.24815]*Cs_scale_factor*Cs_SD_mod;     	
giv_Cs  =   [5.006      3.0724      2.7241      1.3986                  0.6247      0.1972]*Cs_scale_factor;    
U_Cs    =   [];
U_Cs_0  =   []*Cs_scale_factor; 

% gold content in liquid
% M_Cl    =   [Cl2      	Cl3         Cl4         Cl5         Cl6         Cl7         Cl8         Cl9         Cl10        Cl11        Cl12        Cl13        Cl14        Cl16        Cl18        Cl19        Cl21        Cl22];
% M_Cl_0  =   [1.7685   	0.8175      1.8603      2.7946      1.7351      0.55066     3.9791      0.07925     0.4672      1.3097      0.55064     0.2086      0.07923     0.0459      0.07922     0.5506      0.07921     0.0792]*Cl_scale_factor;  
% SD_Cl   =   [0.3537     0.1635      0.37206     0.55892     0.34702     0.11012     0.79582     0.01584     0.09344     0.26194     0.11012     0.04172     0.01584     0.00918     0.01584     0.11012     0.01584     0.01584]*Cl_scale_factor;  
% giv_Cl  =   [1.7283     0.5881      1.834       3.0437      1.6214      0.56431     3.608       0.07935     0.4754      1.3262      0.56434     0.2148      0.07933     0.0456      0.07932     0.5643      0.07931     0.0793]*Cl_scale_factor;  
% U_Cl    =   [Cl15];
% U_Cl_0  =   [1.0328]*Cl_scale_factor;  
% Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
M_Cl    =   [Cl2      	Cl3         Cl4         Cl5         Cl6         Cl7                     Cl9         Cl10                    Cl12        Cl13        Cl14        Cl16        Cl18        Cl19        Cl21        Cl22];
M_Cl_0  =   [1.7685   	0.8175      1.8603      2.7946      1.7351      0.55066                 0.07925     0.4672                  0.55064     0.2086      0.07923     0.0459      0.07922     0.5506      0.07921     0.0792]*Cl_scale_factor;  
SD_Cl   =   [0.3537     0.1635      0.37206     0.55892     0.34702     0.11012                 0.01584     0.09344                 0.11012     0.04172     0.01584     0.00918     0.01584     0.11012     0.01584     0.01584]*Cl_scale_factor*Cl_SD_mod;  
giv_Cl  =   [1.7283     0.5881      1.834       3.0437      1.6214      0.56431                 0.07935     0.4754                  0.56434     0.2148      0.07933     0.0456      0.07932     0.5643      0.07931     0.0793]*Cl_scale_factor;  
U_Cl    =   [Cl15];
U_Cl_0  =   [1.0328]*Cl_scale_factor;

% mass flow of gold
M_Au    =   [];
M_Au_0  =   []*Au_scale_factor;
SD_Au   =   []*Au_scale_factor*Au_SD_mod;
giv_Au  =   []*Au_scale_factor;
U_Au    =   [Au23];
U_Au_0  =   [258.3]*Au_scale_factor;

% Modification 4 & 5 & 6: Properly and smartly Scaling the variables
M_sym   =  	[M_Qs       M_Cw        M_Ql       	M_Cs        M_Cl        M_Au];
M_0     =	[M_Qs_0     M_Cw_0      M_Ql_0      M_Cs_0      M_Cl_0      M_Au_0];
SD      =   [SD_Qs      SD_Cw       SD_Ql       SD_Cs       SD_Cl       SD_Au];
given   =   [giv_Qs     giv_Cw      giv_Ql      giv_Cs      giv_Cl      giv_Au];
U_sym   =  	[U_Qs       U_Cw        U_Ql       	U_Cs        U_Cl        U_Au];
U_0     =	[U_Qs_0     U_Cw_0      U_Ql_0      U_Cs_0      U_Cl_0      U_Au_0];

%% Provided answer
if (test_case_flag == 1)
    given_measured_variables_answer = given;
    length_given_measured_variables_answer = size(given_measured_variables_answer,2);
end

%% Equality constraints function definition
equality_constraints = [Qs1-Qs2;                % Ore mass balance (Balancing the solid only flow rates)
        Qs2-Qs4;
%         Qs4-Qs8;      % Qs8 shall be removed or merged with Qs4   % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs4-Qs10%         Qs4-Qs8;      % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
%         Qs10-Qs11;    % Qs11 shall be removed or merged with Qs10     % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs10-Qs13;%         Qs11-Qs13;  % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs13-Qs16;
        
        % Solution mass balance (Balancing the liquid only flow rates)
        Ql3-Qs2*(1/Cw2-1);
        Qs2*(1/Cw2-1)+Ql7+Ql5-Ql6-Qs4*(1/Cw4-1);
%         Qs4*(1/Cw4-1)-Qs8*(1/Cw8-1);  % Cw8 shall be removed or mergedwith Cw4  % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs4*(1/Cw4-1)+Ql9-Ql5-Qs10*(1/Cw10-1);%         Qs8*(1/Cw8-1)+Ql9-Ql5-Qs10*(1/Cw10-1);  % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
%         Qs10*(1/Cw10-1)-Qs11*(1/Cw11-1);              % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs10*(1/Cw10-1)+Ql14-Ql12-Qs13*(1/Cw13-1);%         Qs11*(1/Cw11-1)+Ql14-Ql12-Qs13*(1/Cw13-1);  % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs13*(1/Cw13-1)+Ql18+Ql17-Ql15-Qs16*(1/Cw16-1);
        Ql12-Ql7-Ql19;
        Ql19+Ql15+Ql20-Ql3;
        Ql6-Ql21;   
        Ql21-Ql22-Ql9-Ql14-Ql18;

        % Physical constraint: Equality of ore concentrations or solid
        % fraction (only for serial streams)
%         Cw8-Cw4;      % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
%         Cw11-Cw10;    % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11

        % Gold mass balance in the ore (solid) and in the solution (liquid)
        Qs1*Cs1+Ql3*Cl3-Qs2*Cs2-Qs2*(1/Cw2-1)*Cl2;      
        Qs2*Cs2+Qs2*(1/Cw2-1)*Cl2+Ql7*Cl7+Ql5*Cl5-Ql6*Cl6-Qs4*Cs4-Qs4*(1/Cw4-1)*Cl4;
%         Qs4*Cs4+Qs4*(1/Cw4-1)*Cl4-Qs8*Cs8-Qs8*(1/Cw8-1)*Cl8;      % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs4*Cs4+Qs4*(1/Cw4-1)*Cl4+Ql9*Cl9-Ql5*Cl5-Qs10*(1/Cw10-1)*Cl10-Qs10*Cs10;%         Qs8*Cs8+Qs8*(1/Cw8-1)*Cl8+Ql9*Cl9-Ql5*Cl5-Qs10*(1/Cw10-1)*Cl10-Qs10*Cs10;    % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
%         Qs10*Cs10+Qs10*(1/Cw10-1)*Cl10-Qs11*Cs11-Qs11*(1/Cw11-1)*Cl11;    % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs10*Cs10+Qs10*(1/Cw10-1)*Cl10+Ql14*Cl14-Ql12*Cl12-Qs13*Cs13-Qs13*(1/Cw13-1)*Cl13;%         Qs11*Cs11+Qs11*(1/Cw11-1)*Cl11+Ql14*Cl14-Ql12*Cl12-Qs13*Cs13-Qs13*(1/Cw13-1)*Cl13;  % Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11, Cs11, Cl11
        Qs13*Cs13+Qs13*(1/Cw13-1)*Cl13+Ql18*Cl18-Ql15*Cl15-Qs16*Cs16-Qs16*(1/Cw16-1)*Cl16;      
        
        % Gold mass balance for the fastest nodes (without leaching)
        Ql12*Cl12-Ql7*Cl7-Ql19*Cl19;
        Ql19*Cl19+Ql15*Cl15-Ql3*Cl3;
        Ql6*Cl6-Ql21*Cl21-Au23;
        Ql21*Cl21-Ql22*Cl22-Ql9*Cl9-Ql14*Cl14-Ql18*Cl18;
    
        % Physical constraint: Equality of gold concentration
        Cl21-Cl22;
        Cl21-Cl9;
        Cl21-Cl14;
        Cl21-Cl18;       
        Cl12-Cl7;
        Cl12-Cl19;
        ];

%% All inequality constraint need to be in =< form. e.g. a+b =< 0
% Note that successive linearization cannot handle inequality constraint
inequality_constraints = [-M_sym';      % Positive value constraint
        -U_sym';

        % Bound constraints
        [Qs1 Qs2 Qs4 Qs8 Qs10 Qs11 Qs13 Qs16]'-100;
        [Ql3 Ql5 Ql6 Ql7 Ql9 Ql12 Ql14 Ql19	Ql21 Ql15  Ql17 Ql18 Ql20 Ql22]'-1000;
        0.2-[Cw2 Cw4 Cw8 Cw10 Cw11 Cw13	Cw16]';
        [Cw2 Cw4 Cw8 Cw10 Cw11 Cw13	Cw16]'-1.0;
        [Cs1 Cs2 Cs4 Cs10 Cs11 Cs13	Cs16 Cs8]'-20;
        [Cl2 Cl3 Cl4 Cl5 Cl6 Cl7 Cl8 Cl9 Cl10 Cl11 Cl12	Cl13 Cl14 Cl16	Cl18 Cl19 Cl21 Cl22 Cl15]'-15;
    
        % Inequality physical constraints
        Cl4-Cl8;
        Cl10-Cl11;
        Cs2-Cs1;
        Cs4-Cs2;
        Cs8-Cs4;
        Cs10-Cs8;
        Cs11-Cs10;
        %Cl10-Cl8;
        Cs13-Cs11;
        Cs16-Cs13];
    
%% Mass balancing equations for calculation of non-redundant variables
non_redundants_equations = [%U_sym(1)+M_sym(2);             % mass balancing equation in form of only equating right part, e.g. equation F1 = F2+F3 shall be written as only F2+F3
                            %M_sym(2)
                            ];  
                        
%% RESULTING COMMENTS ON SOLVING THE GOLD CASE

% AS-it-is basis: 
% Observation: Unobservability of unmeasured variables
% Unmeasured_variables_unobservable =
% [ Cl15, Ql15, Ql17, Ql20]
% NON-OBSERVABILITY: Among them the unobservable variables with independent charateristic are listed below:
% Unmeasured_indpendent_unobservable =
% Ql20
 
% Modification 1: Converting unmeasured variable Ql20 to measured variable
% Observation: Unmeasured_variables_unobservable =
% [ Ql18, Ql17, Ql22]
% NON-OBSERVABILITY: Among them the unobservable variables with independent charateristic are listed below:
% Unmeasured_indpendent_unobservable =
% Ql22

% Modification 2: Converting unmeasured variable Ql22 to measured variable
% Observation: Matrix singularity issue
% the equality constraint values are highly non-scaled
% 1.0e+07 *
% 
%    0.176990000000000
%   -0.000010000000000
%   -0.000010000000000
%   -0.000010000000000
%   -0.000010000000000
%   -0.000010000000000
%   -0.000010000000000
%    0.721941034482759
%    8.342932103516880
%   -0.490008380164880
%    0.003123494981740
%   -0.172260714474428
%   -0.008591262295424
%   -0.767743896112754
%                    0
%   -6.397000000000000
%    0.002000000000000
%   -3.058001000000000
%   -0.000000002640000
%   -0.000000000860000
%    0.000059682086389
%    0.000004284921577
%    0.000000166244773
%    0.000009851116642
%   -0.000000534374930
%   -0.000000356552741
%   -0.000001861615886
%   -0.000000000143500
%   -0.000006953937000
%   -0.000000173481920
%   -0.000000242641489
%    0.000000000000000
%   -0.000000000000000
%   -0.000000000000000
%   -0.000000000000000
%   -0.000000000000000
%    0.000000000000000
% standard_deviation_f_MU = 1.834969272221037e+07

% Modification 3: Properly and smartly Scaling the variables
% We can see that after properly unit conversion the values of Qs and Ql
% are really high. In all the equations, Qs and Ql in all equation are of linear nature. i.e. in each term in all equations:
% 1. Either, Qs or Ql are present
% 2. Qs*Qs or Ql*Ql or Qs*Ql terms that make these terms non-linear are not
% there
% 3. there are no constant terms in these equations
% 4. thus, if these equations which are actually equates to 0 makes no
% difference if altogether they are divided by a constant, and i.e. right
% kind of scaling
% Firstly, only Qs, Ql were scaled that reduced the standard_deviation_f_MU
% to 46.190045213974898

% Modification 4: Au was also scaled as Qs and Ql were scaled, these
% scaling makes the scaling proper linear as in one of the equation Au
% apeared alone with other terms containing Qs or Ql
% observation: standard_deviation_f_MU was reduced to 18.349669701580272
% the trouble causing equations are:
% Qs2*(1/Cw2-1)+Ql7+Ql5-Ql6-Qs4*(1/Cw4-1);
% Ql19+Ql15+Ql20-Ql3;  
% Ql21-Ql22-Ql9-Ql14-Ql18;
% It is probably because of high values of Ql3, Ql6, Ql21

% Cancelled: Modification 5: The max value of Ql i.e. 154.94 is around 2.8 times
% larger than the max value of Qs i.e. 55.8
% Therefore, scaling Ql values by 1/3
% Also, in equation Ql6*Cl6-Ql21*Cl21-Au23; Ql and Au appeared of linear
% nature, thus scaling down Au also by 1/3
% Observation: standard_deviation_f_MU was 20.728868
% Observation: There is a common observation from all of the said changes
% that at the first iteration the SSE value is around 298 while at the
% second it shoots up to 1585 then for further iterations it becomes close
% to 763 then 739 with no further significant changes. At the first
% instance, we know that the SSE has to decrease which is not very apparent
% from the present case. And, therefore, there is some problem. 

% Modification 6: Merging stream 8 into 4 and 11 into 10. Thus, removing
% relevant equations and variables i.e. Qs8, Cw8, Cs8, Cl8, Qs11, Cw11,
% Cs11, Cl11
% Observation: still matrix is close to singularity but now negative value of measured and unmeasured values are not coming




