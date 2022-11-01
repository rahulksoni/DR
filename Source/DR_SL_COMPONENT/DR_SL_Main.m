function [] = DR_SL_Main()
%% This is an auto generated MATLAB file from Optimization Tool.

clear all                   % clear the variables data
close all                   % closes the opened  windows
format SHORTENG
%% Setting up global and local variables 
global M_sym U_sym length_M length_U result input_file_name redundancy_check_status nr_counter nr_positions observability_check_status min_nu_counter equations_repeatibility_check_status variables_nature_repeatibility_check_status equations_repeatibility_flag variables_repeatibility_nature_flag f_MU

% where,
% 1. M_sym is symbolic vector of measured variable
% 2. U_sym is symbolic vector of unmeasured variable
% 3. length_M and length_U is length of M and U vector, respectively
% 4. result is the final matrix to save row wise results in specified sequence 
% 5. input_file_name is the name of the matlab .m file that holds all the case
% to case input data
% 6. redundancy_check_status is flag that says the case is checked for
% presence of redundant variables, this is required because redundancy
% check has to be performed only once
% 7. nr_counter shall count the number of non-redundant variables
% 8. nr_positions returns the position in M_sys of non-redundant variables
% 9. observability_check_status is flag that says the case is checked for
% the presence of unobservable unmeasured variables, this is required because
% observability check has to be performed only once
% 10. min_nu_counter shall count the minimum number of non-observable
% variables, why we are calling minimum here. Because as per theory, the
% non-zero rows in inv(R1)*R2 matrix decides whether the calculation of rest unmeasured
% variables depend upon these minimum number of unobservables. If depends,
% then definetely the total number of unobservables will increase.
% 11. equations_variables_repeatibility_check_status is the trigger that whether equality
% constraints have been checked for repeatitions in equations
% 12. equations_repeatibility_flag turns from 0 to 1 if there are repeated
% equality constraints
% 13. variables_repeatibility_nature_flag turns from 0 to 1 if there are
% linearly dependent columns in combined incidence matrix
% 14 f_MU stores the value vector of values of equality constraints

%% Input data/material file content/text need to be copied in the file inputs.m
% The content/text of file name for input data updated in file inputs.m is imported in program through command eval
input_file_name = 'DR_SL_Inputs';
eval(input_file_name);   

a = sl;
if (a == sqp)
    run DR_Main_SQP.m
elseif (a == sl)
    run run DR_Main_SL.m
end

%% Setting up font sizes that are going to be used in plots
titlefontsize = 18;
axisfontsize  = 16;
labelfontsize = 18;

%% Global variable settings  
length_M = size(M_sym,2);                           % number of measured variables
length_U = size(U_sym,2);                           % number of unmeasured variables

%% Some calculations
V = power(SD,2);                                    % Variance vector for the measured variables
W = diag(V);                                        % Co-variance matrix 

%% Redundancy (measured variables) and observability (unmeasured variables) related
redundancy_check_status = 0;                        % initially it is 0; once redundancy check is done it will turn to 1
nr_counter              = 0;                        % counter to count the number of non-redundant variables
observability_check_status = 0;                     % initially it is 0; once observability check is done it will turn to 1
min_nu_counter          = 0;                        % counter to count the number of non-observable variables
equations_variables_repeatibility_check_status = 0; % just a trigger that ensure once to be performed repeatibility check of equality constraints is done 
equations_repeatibility_flag = 0;                   % turns to 1 if there are linearily dependent equality constraints
variables_repeatibility_nature_flag = 0;            % turns to 1 if there are linearily dependent linearily dependent columns in combined incidence matrix

%% Some pathological checkups on the data provided
if ((test_case_flag == 1) && (length_given_measured_variables_answer ~= length_M)) 
    disp(['Error: the number of measured variables in given answer is not equal to number of measured variables'])
    return
end
if (size(M_0,2) ~= length_M) 
    disp(['Error: the initials for measured variables are not equal to number of measured variables'])
    return
end 
if (size(U_0,2) ~= length_U) 
    disp(['Error: the initials for unmeasured variables are not equal to number of unmeasured variables'])
    return
end 
if (size(V,2)   ~= length_M) 
    disp(['Error: the SD for measured variables are not equal to number of measured variables'])
    return
end 

%% Checking for the exact duplicate values in M_0 and U_0
% Before that, why U_0 is required. Unmeasured variables need to be
% predicted then why its initial values is required?
% Answer, is because we are going to check if there are dependent columns
% or variables in the incidence matrix A or [Am Au]. After, taking Jacobian
% of bilinear system there are many variables that get removed from the 
% column or row making them linear column or row. In such case, unmeasured
% variables also play the role and their duplication need to be checked.
% Secondly, if we see the following part in ObjectiveFunction.m
% Am      = double(subs(Am_jacobian));                        
% Au      = double(subs(Au_jacobian)); 
% Then, we can understand that why real values of the unmeasured variables
% at the initial stage is required. 
% There might be linearily dependent
% columns in incidence matrix or A or [Am Au]. There are various possible
% reasons for the same. One for the obvious reason comes from the fact
% that in successive linearization scheme of bilinear equality constraint,
% the equations are partially derivated by each variable through jacobian
% scheme. Partial derivation removes many variables from equation resulting
% in shorter equation or sometimes single variable answer. Now, if there
% are different variables of same initial value then there are chances that
% that due to same values columns become linearily dependent. Therefore, the
% initial values in M_0 and U_0 combinely must be all unique (even by slight difference).
no_duplicate_found_flag = 0;    % shall turn to 1 if there are no duplicates found
while no_duplicate_found_flag < 1
    duplicate_initials_flag = 0;    % shall turn to 1 if there are exact duplicate initial values
    M_0_mod = M_0;                  % temporary allocation
    U_0_mod = U_0;                  % temporary allocation
    len_max = length(M_0_mod);      % temporary allocation
    if (length_M > length_U)
        U_0_mod = [U_0_mod NaN(1,length_M-length_U)];   % matching the length of U_0 with M_0 with NaN values
    else
        M_0_mod = [M_0_mod NaN(1,length_U-length_M)];   % matching the length of M_0 with U_0 with NaN values
    end
    len_max  = length(M_0_mod);     % now both measured and unmeasured have same length that equats to len_max
    M_0_mod2 = num2cell(M_0_mod);   % allocating as cell values so that duplicate positions can be replace with character 'D'
    U_0_mod2 = num2cell(U_0_mod);   % allocating as cell values so that duplicate positions can be replace with character 'D'
    SD_mod   = num2cell(SD);        % allocating as cell values so that duplicate positions can be replace with character 'D'
    if (M_0_mod(1) == U_0_mod(1))   % comparing the first values of M_0_mod and U_0_mod. following loop shall not compare first values and that is why this is required
        M_0_mod2{1} = 'D';          % if matches then replace first of M_0_mod2 with 'D'
        % we can slightly modify the duplicate number by additing small
        % number but we have to do several times for other numbers so the
        % number to be added has to be randomly generated number. rand
        % command in matlab generates random numbers between 0 and 1 which
        % is wide range. The random numbers can be generated between 1
        % and 5 by (1+(5-1)*rand(1)) which has to be scaled down to
        % small value may be by the order of 100 to give % value change.
        M_0(1)      = M_0(1)*(1+(1+(5-1)*rand(1))/100); % making the value unique by adding small random number 
        duplicate_initials_flag = 1;    
    end
    for i=1:len_max       
        for j=i+1:len_max
            if ((M_0_mod(i) == M_0_mod(j)) || (M_0_mod(i) == U_0_mod(j)))   % comparing measured variable 1 by 1 to all other measured and unmeasured variables and marking them with 'D' any duplicate of the same is found
                M_0_mod2{i} = 'D';
                M_0(i)      = M_0(i)*(1+(1+(5-1)*rand(1))/100);         % making the value unique by adding small random number 
                duplicate_initials_flag = 1;
            end
            if ((U_0_mod(i) == U_0_mod(j)) || (U_0_mod(i) == M_0_mod(j)))   % now, comparing unmeasured variable 1 by 1 to all other unmeasured and measured variables and marking them with 'D' any duplicate of the same is found
                U_0_mod2{i} = 'D';
                U_0(i)      = U_0(i)*(1+(1+(5-1)*rand(1))/100);         % making the value unique by adding small random number 
                duplicate_initials_flag = 1;
            end
        end
    end
    if (duplicate_initials_flag == 1)
        S1 = sprintf('%.1s  ', M_0_mod2{:});
        S2 = sprintf('%.1s  ', U_0_mod2{:});
        disp(['Error-DUPLICATE-INITIALS: there are initial values in the initial value of measured (M_0) and unmeasured (U_0) variables that exactly matches other than its own value'])
        disp(['M_0 values marked with D for duplicate values'])
        disp(S1)
        disp(['U_0 values marked with D for duplicate values'])
        disp(S2)
        disp(['Error-DUPLICATE-INITIALS: attempts have been made to slightly modify them to make them unique'])
        % return
    end
    if (duplicate_initials_flag == 0)
        no_duplicate_found_flag = 1;
    end
end

%% Recheck to see if the values have been modified to make them unique
duplicate_initials_flag = 0;
disp(['Duplicate values have been slightly modified to make them unique'])

%% Zero value of standard deviation create problems as observed in one of the case. Checking to see if there are 0 values of standard deviations
SD_0_flag = 0;                                      % shall turn to 1 if there are standard deviation values of 0 for the measured variables
SD_0_positions = [];                                % a nill vector to be filled later
for i=1:length(SD)
    if (SD(i) == 0)
        SD_0_positions(size(1,2)+1) = i;
        SD_0_flag = 1;
    end
end
if (SD_0_flag == 1)
    disp('Error: there are standard deviation values of 0 at following position of measured variables, give a smaller for safe side')
    SD_0_positions
    return
end
    
%% Creating result matrix and filling it up with initial values
result(1,1:length_M)                    = M_0;      % initial values of measured variables
result(1,length_M+1:length_M+length_U)  = U_0;      % initial guessed values of unmeasured variables
result(1,length_M+length_U+1)           = 0;        % initial Obj_fun
result(1,length_M+length_U+1+1)         = 0;        % absolute of difference between consequtive measured variables (y) for convergence 
result(1,length_M+length_U+1+1+1)       = 0;        % absolute of difference between consequtive unmeasured variables (z) for convergence 
result(1,length_M+length_U+1+1+1+1)     = 0;        % average value of equality contraint function or fun_MU or f_MU; ideally it has to be zero. So, if this is reducing then we are in the right direction
result(1,length_M+length_U+1+1+1+1+1)   = 0;        % average percent value of absolte of deviations from measured values of measured variables 

%% External function calls, their jacobian expressions
fun_MU      = equality_constraints;                 % symbolic equality constraints, the complete incidence matrix for finding the jacobian as Am and Au
Am_jacobian = jacobian(fun_MU,M_sym);               % symbolic incidence matrix for measured variables
Au_jacobian = jacobian(fun_MU,U_sym);               % symbolic incidence matrix for unmeasured variables
m           = size(Au_jacobian,1);                  % number of rows in Au or say the number of equality constraints
n           = size(Au_jacobian,2);                  % number of columns in Au or say the number of unmeasured variables

%% number of columns in the Au_jacobian has to match no of unmeasured variables    
% if (m <= n)                                          % always, m > n
%    disp(['Error: the number of independent mass balancing equations (equality constraints) are not sufficient; dropping off calculations'])
%    disp(['Explanation: As per theory, if Au (separated incidence matrix for unmeasured variables) is an mxn matrix then in qr factorization, Q and R has to be mxm and mxn matrix, respectively. However, some of the rows in R are always 0, eventually resulting in R1 as nxn top positioned upper triangular matrix of R.']) 
%    disp(['Q1 and Q2 are the left and right partitions of the Q matrix, respectively. Here, Q1 is left mxn matrix and Q2 is the remaining right mx(m-n) matrix.'])
%    disp(['Simple, implication is that for the existence of m > n'])
%    disp(['Here, m and n are also equal to number of equality constraints and number of unmeasured variables, respectively'])
%    return
% else
if (length_U ~= n)
   disp(['Error: number of columns in partitioned incidence matrix for unmeasured variables is not equal to number of unmeasured variables; dropping off calculations'])
   return
end

%% Following section is not required as Optimization command is not used in successive linearization scheme
% %% Start with the default options
% options = optimoptions('fmincon');
% %% Modify options setting
% options = optimoptions(options, 'MaxFunEvals', 10000, 'TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-6);
% options = optimoptions(options,'Display', 'iter-detailed');
% options = optimoptions(options,'Algorithm', 'sqp');
% [x,fval,exitflag,output,lambda,grad,hessian] = ...
% fmincon(@LagrangeFunction,[y_0],[],[],[],[],[],[],@ConstraintFunction,options);

%% Combining both measured unmeasured symbolic variable names for use in printing variables names in resulting tables 
% such exercise is also required as symbolic names cannot be used for xtick mark in one of the plot
sym_variable_names = [M_sym U_sym];
sym_variable_names_transposed = sym_variable_names.';
for i=1:(length_M+length_U)  
    temp1 = sym_variable_names_transposed(i);                   % getting the variable's symbolic name
    temp2 = char(temp1);                                        % converting name to character form
    x_axis_marks_measured_unmeasured(i,1) = {temp2};            % these tick mark names shall be used in one of the iterative plot as xaxis names
end
M_sym_transposed = M_sym.';
for i=1:length_M 
    temp1 = M_sym_transposed(i);                                % getting the variable's symbolic name
    temp2 = char(temp1);                                        % converting name to character form
    x_axis_marks_measured(i,1) = {temp2};                       % these tick mark names shall be used in one of the iterative plot as xaxis names
end

%% Stats printing
fprintf('Problems stats:')
fprintf('\nNumber of equality constraints = %d',m)
fprintf('\nNumber of measured variables   = %d',length_M)
fprintf('\nNumber of unmeasured variables = %d\n',length_U)

%% Copyright preventer
crs = DR_SL_crsc();if (crs == 1),    return,    end

%% Printing outputs to a file
% Cleaning the text inside file
f_name = fopen('outputParameters&Results.csv','w' );
    fprintf(''); 
fclose(f_name);
% append case output data in cleaned file
f_name = fopen('outputParameters&Results.csv','a+' );
fprintf(f_name, 'Problems stats:');
fprintf(f_name, '\nNumber of equality constraints = %d',m);
fprintf(f_name, '\nNumber of measured variables   = %d',length_M);
fprintf(f_name, '\nNumber of unmeasured variables = %d\n',length_U);

%% while loop shall be used until the convergence criteria is fulfilled 
figure('units','normalized','outerposition',[0 0 1 1])
hold on
while convergence_flag < 1
    previous_redundancy_check_status = redundancy_check_status;                         % recording the redundancy check status before functions runs
    [SAE_M SAE_U] = DR_SL_ObjectiveFunction(W,fun_MU,Am_jacobian,Au_jacobian,m,n,non_redundants_equations,zero_tol,x_axis_marks_measured_unmeasured,f_name,titlefontsize,axisfontsize,labelfontsize);            % sum of absolute error for last two iteration
    if (equations_repeatibility_flag == 1)
        disp(['Error-EQUATIONS_REPEATIBILITY: there are linearily dependent equality constraint; consider removing them to avoid matrix singularity error'])
        break
    end    
    if (min_nu_counter > 0)
        disp(['Error-NON_OBSERVABILITY: consider removing unobservable unmeasured variables ONE-BY-ONE; starting from the unobservables variable with independent characteristic, either by converting them to measured variables or by providing additional mass balancing euqations'])
        disp(['Error-NON_OBSERVABILITY: such removal is necessary because, in successive linearization scheme, the new values of unmeasured values are also used for the estimation of jacobian in the next iteration and therefore their numeric values are required'])
        break
    end 
    if (redundancy_check_status - previous_redundancy_check_status == 1)        % ensured that the redundanct check is performed just now and it is passed; stopping here to ask for manual inputs on reconciliation of non-redundant variables       
        if ((nr_counter > 0) && (length(non_redundants_equations) ~= nr_counter))   % checking if number of non-redundant measured variables is more than 0 and sufficient number of equation are not provided to recalculate them
            disp(['CORRECTION-NON_REDUNDANCY: NON-REDUNDANT variables have been found, please provide appropriate mass balancing equations in input file for their recalculation.'])
            disp(['CORRECTION-NON_REDUNDANCY: The additional equations MUST NOT INCLUDE NON-OBSERVABLE unmeasured variables.'])
            disp(['CORRECTION-NON_REDUNDANCY: The mass balancing equation must be in the trasposed form of above order.'])
            disp(['Error-NON_REDUNDANCY: The number of mass balancing equations provided for recalculating non-redundant variables are not equal to the number of non-redundant variables'])
            break
        elseif (((nr_counter == 0) && (length(non_redundants_equations) > 0)))
            disp(['Error-NON_REDUNDANCY: Even though there are no non-redundant variables; the non-redundant equations have been provided. Consider removing or commenting the equations.'])
            break
        end        
    end
    if ((SAE_M < convergence_tol_M) && (SAE_U < convergence_tol_U))
        convergence_flag = 1;
    end
end 
hold off

%% if, there were issues in the while that broke it, then there is no point in going ahead
if (convergence_flag == 0)
    disp(['Errors-occurred: terminating the program'])
    return
end
 
%% Results section
n_rows_result = size(result,1);

%% following 2 loops are required to do two things: first to convert symbolic variables names to character variables names and secondly all character varaibles names of same size
%% otherwise, the character variables of different size shall create problem in next step of printing variable names and numeric answer from result matrix
max_char_size = 0;
for k=1:(length_M+length_U)
    temp1   = sym_variable_names_transposed(k);        % temporary allocation of a symbolic variable
    temp2   = char(temp1);              % convertingt allocated variable to character variables
    temp3   = length(temp2);            % getting the length of it
    if (max_char_size < temp3)          % loop to get max size of biggest variable name
        max_char_size = temp3;          
    end 
end

% equalizing the size of variable names
extra_blank_space = 5;                  % extra blank space shall be appended at the end of the character variable or later convenient look in matrices
for k=1:(length_M+length_U)
    temp1   = sym_variable_names_transposed(k);    % temporary allocation of a symbolic variable
    temp2   = char(temp1);              % convertingt allocated variable to character variables
    temp3   = length(temp2);            % getting the length of it
    b       = blanks(max_char_size - temp3 + extra_blank_space);    % creating blank space of difference of max to current size + extra blank space for convenient look at later stage
    V(k,1:(max_char_size + extra_blank_space)) = [temp2 b];  % concatenating actual variable with blank spaces
    V       = char(V);                  % converting the variable from coded to character form
end

% equalizing the size of term 'Names' to previous common size of variable
% names
temp1   = char('Names');                % shall be one of the title of table
temp2   = length(temp1);                % length of term Variable_names
b       = blanks(max_char_size - temp2 + extra_blank_space);    % blanks required to convert Variable_names to the size of previous V
Names   = [temp1 b];                    % Names with added spaces
Names   = char(Names);                      

%% Variables iterative progress table printing
result_measured_transpose       = result(:,1:length_M)';
result_unmeasured_transpose     = result(:,length_M+1:length_M+length_U)';
result_measured_unmeasured      = [result_measured_transpose;result_unmeasured_transpose];
n_iterations                    = size(result_measured_transpose,2);        % total number of iteration the program was ran
iteration_vector                = 1:n_iterations;

%% Printing the final results to matlab command prompt and to a file
% printing
fprintf('\nVariables iterative progress')
fprintf(f_name, '\nVariables iterative progress');
fprintf('\n%s',Names)
fprintf(f_name, '\n%s',Names);
for k=1:n_iterations
    fprintf('iter%d         ',iteration_vector(k))
    fprintf(f_name, 'iter%d         ',iteration_vector(k));
end
for k = 1:(length_M+length_U)
    variable_value_vector = sprintf('%f      ', result_measured_unmeasured(k,:));
    fprintf('\n%s%s     ',V(k,:),variable_value_vector)
    fprintf(f_name, '\n%s%s     ',V(k,:),variable_value_vector);
end
 
%% Variables and deviations final values table printing 
% initial values
Measured_initial_values              = M_0;
Unmeasured_initial_values_fake       = NaN(1,length_U);                % this won't contain the U_0 initial values as those are not provided by plant as for measured variables, they are just starting point particularly required by this  method
Measured_Unmeasured_initial          = [Measured_initial_values Unmeasured_initial_values_fake];

% obatined values
Measured_variables_obtained_values   = result(n_rows_result,1:length_M); 
Unmeasured_variables_obtained_values = result(n_rows_result,length_M+1:length_M+length_U);
Measured_Unmeasured_obtained         = [Measured_variables_obtained_values Unmeasured_variables_obtained_values];   % containing answers for both measured and unmeasured variables

% deviations
deviation_percent_measured_variables = (Measured_variables_obtained_values - M_0)./M_0*100;           % calculates the percent deviation in obtained measured variables from M_0
deviation_percent_unmeasured_fake    = NaN(1,length_U);
deviations_percent                   = [deviation_percent_measured_variables deviation_percent_unmeasured_fake];

% printing 
fprintf('\n\nFinal results with deviations')
fprintf(f_name, '\n\nFinal results with deviations');
fprintf('\n%sinitial       obtained      deviations(%%)\n',Names)
fprintf(f_name, '\n%sinitial       obtained      deviations(%%)\n',Names);
for k=1:(length_M+length_U)
    fprintf('%s%8f      %8f      %8f\n',V(k,:),Measured_Unmeasured_initial(k),Measured_Unmeasured_obtained(k),deviations_percent(k))    % number Between % and f decides the length of printed variable, it is used becuase real values and NaN have different size 
    fprintf(f_name, '%s%8f      %8f      %8f\n',V(k,:),Measured_Unmeasured_initial(k),Measured_Unmeasured_obtained(k),deviations_percent(k));    % number Between % and f decides the length of printed variable, it is used becuase real values and NaN have different size 
end

%% Obj_fun values printing
Obj_fun_obtained_value = sum(power((Measured_variables_obtained_values - M_0),2));
if (test_case_flag == 1)    
    Provided_measured_variables_Obj_fun_value = sum(power((given_measured_variables_answer - M_0),2));   
    fprintf('\n\nProvided & obtained measured variable Obj_fun values:')
    fprintf(f_name, '\n\nProvided & obtained measured variable Obj_fun values:');
    fprintf('\n%8f   %8f\n',Provided_measured_variables_Obj_fun_value,Obj_fun_obtained_value)
    fprintf(f_name, '\n%8f   %8f\n',Provided_measured_variables_Obj_fun_value,Obj_fun_obtained_value);
end

%% Plotting the final resulting values and their deviations
if (test_case_flag == 1)   
    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot (2,1,1)
    hold on 
    %plot(given_measured_variables_answer)
    plot(given_measured_variables_answer,'ro',...
        'LineWidth',2,...
        'MarkerSize',4,...
        'MarkerEdgeColor','red')    
    plot(Measured_variables_obtained_values,'s',...
        'LineWidth',2,...
        'MarkerSize',8,...
        'MarkerEdgeColor','green')  
    set(gca,'xtick',1:length_M)
    set(gca,'xticklabel',x_axis_marks_measured)
    set(gca,'XTickLabelRotation',45)
    legend('provided values','obtained values') 
    legend('Location','best')
    title('Provided and obtained values for measured variables','FontSize',titlefontsize)
    set(gca,'FontSize',axisfontsize)
    ylabel('values','FontSize',labelfontsize)
    grid on
    grid MINOR     
    hold off
    
    subplot (2,1,2)
    hold on
    line_45 = linspace(min(given_measured_variables_answer),max(given_measured_variables_answer),100);
    plot(line_45,line_45,'--')
    plot(given_measured_variables_answer,Measured_variables_obtained_values,'ro',...
        'LineWidth',2,...
        'MarkerSize',8)
    legend('ideal fall-off','obtained fall-off')
    legend('Location','southeast')
    title('Correlation plot for measured varaibles','FontSize',titlefontsize)
    set(gca,'FontSize',axisfontsize)
    xlabel('provided values','FontSize',labelfontsize)
    ylabel('obtained values','FontSize',labelfontsize)
    grid on
    grid MINOR
    hold off
end

%% Plotting the deviation of measured variables
% errorbar plot, plots the provided error value on both sides i.e. up and down
% so, if half the deviations then these half magnitude shall be plotted on up and down side of this mid value
% this shall fulfill our ultimate objective of plotting total deviation or
% say actual unhalved deviation
zero_vector     = zeros(length_M);                                  % to create a horizontal plot at 0 for deviation to be seen in + and - direction
max_deviation   = max(abs(deviation_percent_measured_variables));   % max deviation is required to set the y-axis min and max range
half_deviations = (deviation_percent_measured_variables)/2;         % as per descision total deviation shall be plotted in up and down equally of half deviation
y_max           = 1.2*max_deviation;                                % 20% absolute increase in plot y-axis min and max values
x               = 1:length_M;                                       % shall be x-axis hidden values
err             = abs(half_deviations);                             % error values as required by error bar plot. this are absolute of half of the error absolute deviations

figure('units','normalized','outerposition',[0 0 1 1])              
hold on
plot(zero_vector,'--g','LineWidth',1)                               % plot the horizontal zero line
errorbar(x,half_deviations,err,'b.','LineWidth',3)                  % error bar plot with given set of conditions
set(gca,'xtick',1:length_M)                                         % setting up hidden x-values
set(gca,'xticklabel',x_axis_marks_measured)                         % setting up visible mask values in x-axis
set(gca,'XTickLabelRotation',45)                                    % rotating the mask value in 45b degrees to avoid overlapping
ylim([-y_max y_max])                                                % setting up min max values of y-axis in plot
title('Percent deviation of measured variables','FontSize',titlefontsize) 
grid on
grid MINOR
set(gca,'FontSize',axisfontsize)
ylabel('deviations','FontSize',labelfontsize)
hold off

fclose(f_name);
f_MU
%% Writing generated data to the file name xlfilename with specified sheet name 
x_write        = xlswrite(xlfilename,Measured_Unmeasured_obtained',SLSheetName,x_writerange)';
equality_write = xlswrite(xlfilename,f_MU,SLSheetName,equality_write_range)';
final_values = Measured_Unmeasured_obtained'
equality_function_values = f_MU