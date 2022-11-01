function [SAE_M SAE_U] = ObjectiveFunction(W,fun_MU,Am_jacobian,Au_jacobian,m,n,non_redundants_equations,zero_tol,x_axis_marks_measured_unmeasured,f_name,titlefontsize,axisfontsize,labelfontsize)
%% Summary of this function goes here

%% Setting up global and local variables 
global M_sym U_sym length_M length_U result input_file_name redundancy_check_status nr_counter nr_positions observability_check_status min_nu_counter equations_repeatibility_check_status variables_nature_repeatibility_check_status equations_repeatibility_flag variables_repeatibility_nature_flag f_MU
M_sym = M_sym; U_sym = U_sym; length_M = length_M; length_U = length_U; result = result; redundancy_check_status = redundancy_check_status; nr_counter = nr_counter; nr_positions = nr_positions; observability_check_status = observability_check_status; min_nu_counter = min_nu_counter; equations_repeatibility_check_status = equations_repeatibility_check_status; variables_nature_repeatibility_check_status = variables_nature_repeatibility_check_status; equations_repeatibility_flag = equations_repeatibility_flag; variables_repeatibility_nature_flag = variables_repeatibility_nature_flag;

%% Getting the values from last iteration based updates result matrix
M_0             = result(1,1:length_M);
n_rows_result   = size(result,1);                                           % the counting till last iteration
M               = result(n_rows_result,1:length_M);                         % values of measured variables from last iteration 
U               = result(n_rows_result,length_M+1:length_M+length_U);       % values of unmeasured variables from last iteration

%% Filling up symbolic variables with previous iteration numeric values for subsitutions in jacobian matrices at a later stage
for i = 1:length_M
    evalc(sprintf('%s=M(%d)',char(M_sym(i)),i));            
end
for i = 1:length_U
    evalc(sprintf('%s=U(%d)',char(U_sym(i)),i));
end

%% Function real value substitutions
f_MU    = double(subs(fun_MU));                             % the equality constraints with real values
Am      = double(subs(Am_jacobian));                        % incidence matrix for measured variables with real values
Au      = double(subs(Au_jacobian));                        % incidence matrix for unmeasured variables with real values

%% Getting the full incidence matrix for checkup purpose only
A   = [Am Au];                                              % the complete linearized form of incidence matrix

%% Lets check if there are repeating mass balance equations (equality constraints)
    % the linearily dependent equality constraints shall create the matrix
    % singularity issue and therefore need be removed. sicne, each row in 
    % the combined incidence matrix represent corresponding equality
    % constraint, the linearily dependent rows and correspongding 
    % equations shall required to be removed
if (equations_repeatibility_check_status == 0)
    [n_ind_row_A,ind_col_vector_AT,AT_reduced] = linearColumnDependency(A');        % The program linearColumnDependency checks for linear dependent columns then produces position of independent columns and also the reduced matrix after eliminating dependent columns. AT means A transposed.
    ind_row_vector_A = ind_col_vector_AT';
    if (n_ind_row_A < m)                                                            % where, m is the number of equality constraints provided
        n_linearily_dependent_equations = m - n_ind_row_A;
        
        % Now, we will list the dependent equation positions, in a two-step
        % process; first listing all the position and then removing
        % independent ones
        index_dependent_equality_constraint                     = (1:m)';           % making a vector of position of all equality constraints
        index_dependent_equality_constraint([ind_row_vector_A]) = [];               % then, removing the independent ones from them
             
        fprintf('\nEQUATIONS-REPEATIBILITY: Consider removing the %d linearily dependent equality constraints out of total %d equality contraint.\n',n_linearily_dependent_equations,m)
        disp(['EQUATIONS-REPEATIBILITY: the positions of dependent equality constraints  are as follows (after removing one independent instance from each group of dependent equations):'])
        fprintf('%d\n',index_dependent_equality_constraint)   

        equations_repeatibility_flag = 1;
        SAE_M = 0; SAE_U = 0;   % temporary values just to finish loop without errors
        return
    end
    
    equations_repeatibility_check_status = 1;
end

%% important points about the rank of matrix
% 1. A square matrix is a full rank if its determinant is non-zero (If determinant is zero or close to zero then the matrix is not full rank)
% 2. rank deficiency in matrix is caused by row-row or column-column linear dependence 
% 3. In mxm matrix, the row rank and column rank are always equal
% 4. In mxn matrix, max possible rank is m (if, m < n) or n (if, n < m)
% 5. In non-square matrix the max possible rank is number of row, if its
% number is leeser than number of columns and vice versa
% 6. The rank of a matrix directly relates to the number of independent
% vectors in the direction of lowest length/dimension
%% Some checkups regarding unobservability of unmeasured variables
% As per the theory, if Au is an mxn matrix then in qr factorization, Q and R has to be mxm and nxn
% matrix (with other bottom rows in R as zero). As such, matlab will give R
% in mxn size but since bottom rows will be zero, it is eventually as nxn matrix
% The formulation of Q1 and Q2 (Ref: http://www.polymtl.ca/namp/docweb/Modules_Web/M11_Tier1_Chap1-3.pdf)
% As per theory the Q1 is mxn sized left part of Q matrix and remaining
% that is right mx(m-n) part is Q2. And, R1 is the nxn top positioned upper triangular matrix of R.
if (observability_check_status == 0)
    if (rank(Au) == n)        
        fprintf('\nNON-OBSERVABILITY: All the unmeasured variables are OBSERVABLE\n')
    elseif (rank(Au) < n)
        fprintf('\nNON-OBSERVABILITY: the partitioned incidence matrix for unmeasured variables is not a full rank matrix\n')
    end
    [n_ind_col,ind_col_vector,Au_reduced] = DR_SL_linearColumnDependency(Au);         % The program linearColumnDependency checks for linear dependent columns then produces position of independent columns and also the reduced matrix after eliminating dependent columns 
    if (n_ind_col < n)
       min_nu_counter = n - n_ind_col;                                          % the minimum number of unobservables unmeasured variables; minimum because other unmeasured variables may depend upon these, as checked by inv(R1)*R2 matrix
       disp(['NON-OBSERVABILITY: there are MINIMUM ', num2str(min_nu_counter), ' out of ', num2str(n), ' unmeasured variables that are unobservable'])  
    end    
end

%% Most imp n value that is used for forming Q1, Q2 and R1, R2 need to be changed if there are unobservable unmeasured variables
       % Now as per theory the size of column wise partition of Q in Q1 and
       % Q2 need to be updated by effectively reducing the column size of
       % Q1. For example, if there 35 (n=35) unmeasured variables and then there
       % are only 24 independent columns in the Au. Then the previous value
       % of 31 shall change to 24. There were 31 equality contraints. Means
       % Q size is 31x31, now how come we Q1 of size 31x35 out of 31x31 Q.
       % Since, there are 24 indpenedent column in Au, the effective n=24.
       % Now, one can split the Q (31x31) in Q1 (31x24) & Q2 (31x(31-24)).
       % For full details see file M11_Tier2_case1.pdf
if (min_nu_counter > 0)                                                         % if non-observable unmeasured variables exist
    n_reduced = n - min_nu_counter;
end

%% Primary computation code
b           = Am*M(1,:)'+Au*U(1,:)'-f_MU;

if (min_nu_counter == 0)                        % if there are no dependent columns in Au or say no unobservable unmeasured variables. Formation of Q1, Q2, R1, R2 depends upon the rank of Au.              
    [Q,R,E] = qr(Au);
    Q1      = Q(1:m,1:n);
    Q2      = Q(1:m,n+1:m);
    R1      = R(1:n,1:n);
    inv_R1  = inv(R1);
else                                            % if there exist unobservable unmeasured variables
    [Q,R,E] = qr(Au);
    Q1      = Q(1:m,1:n_reduced);
    Q2      = Q(1:m,n_reduced+1:m);
    R1      = R(1:n_reduced,1:n_reduced);
    R2      = R(1:n_reduced,n_reduced+1:n);
    inv_R1  = inv(R1);
    Ru      = inv(R1)*R2;
    
    %% following are different property based positions of U_sym vector
    U_positions_original            = [1:1:length_U]';
    U_positions_permuted            = E'*U_positions_original;                                       % refer to theory in notebook that permuted u is equal to E'*u
    U_n_positions                   = U_positions_permuted(1:end-min_nu_counter);                    % these positions are of unmeasured variables of which some may possibly be unobservable 
    U_positions_n_n_reduced         = U_positions_permuted(end-min_nu_counter+1:end);                % these positions of unmeasured variables are surely unobservable as counted by min_nu_counter
    % Now, we will see which rows in the Ru is all 0. find(all(abs(Ru) <
    % zero_tol, 2)) command returns the position of rows in Ru that are all
    % 0. Then, we will take those positions in Ru and variables at same
    % position in U_sym shall become U_n_observable_positions, as these variables do
    % not depend on already unobservable variables as dictated by positions
    % in U_positions_n_n_reduced of U_sym variables.
    % Rest variables dictated by U_n_unobservable_positions in U_sym shall be unobservable.                                                 % when these are no dependeny i.e. (find(all(abs(Ru) < zero_tol, 2))) return no value then next line throws some problem. to safeguard this if else condition is designed
    U_n_observable_positions        = U_n_positions(find(all(abs(Ru) < zero_tol, 2)));               % these are unmeasured variables that are concurrently unobservable due to dependency
    U_n_unobservable_positions      = U_n_positions(find(all(abs(Ru) > zero_tol, 2)));               % results the positions of from U_n except U_n_observable positions; these are non-observable concurrents non-observable positions depending each-other
    U_sym_unobservable_positions    = [U_n_unobservable_positions; U_positions_n_n_reduced];         % the U_sym positions of all unobservables variables
    U_sym_observable_positions      = U_n_observable_positions;                                      % the U_sym positions of observale variables 
end

%% Now one time printing of non-observable variables for information/improvement purpose
if (observability_check_status == 0)
    if (min_nu_counter > 0)
        U_observables_unobservable = U_sym;                                      % trial allocation for modification purpose at later stage
        U_observables_unobservable(U_sym_unobservable_positions) = 'UNOBSERVABLE'; % shall print UNOBSERVABLE in place of unobservable variables
        fprintf('NON-OBSERVABILITY: The original/provided unmeasured variable vector is:\n')
        Unmeasured_variables_original       = U_sym                             % just to show original unmeasured variables in given sequence
        disp(['NON-OBSERVABILITY: Among them the unobservable variables are marked below:'])
        Unmeasured_variables_observability  = U_observables_unobservable        % shows the unobservable variable by marking them as 'UNOBSERVABLE'
        disp(['NON-OBSERVABILITY: They are collectively listed below:'])
        Unmeasured_variables_unobservable   = U_sym(U_sym_unobservable_positions)        % shows the unobservable variables
        disp(['NON-OBSERVABILITY: Among them the unobservable variables with independent charateristic are listed below:'])
        Unmeasured_indpendent_unobservable  = U_sym(U_positions_n_n_reduced)    % shows only the independent unobservable unmeasured variables, that are actually made by U_positions_n_n_reduced
        SAE_M = 0; SAE_U = 0;                                                   % allocating trial values just to end the loop without error
        return
    end
    observability_check_status = 1;                                             % this says that the on time redundancy check is performed
end

%% till this point we were only dealing with the jacobian or partial derivative of the equality constraint function, now we are going to actually deal with the derivated equation that arised from the condition that minimize[(M-Mhat)*'inv*(V)(M-Mhat)] when constrained by Am*M+Au*U=b
%% Calculation for new values of measured variables
Q2TAmWQ2TAmT        = (Q2' * Am) * W * (Q2' * Am)';
inv_Q2TAmWQ2TAmT    = inv(Q2TAmWQ2TAmT);                                                                % getting inverses at hand for error debugging
right_second_part   = W * (Q2' * Am)' * inv_Q2TAmWQ2TAmT * (Q2' * Am * M(1,:)' - Q2' * b);              % separate calculation of this part is necessary to know if measured variables are non-redundant
Mhat                = M(1,:)' - right_second_part;                                                      % getting reconciled measured variables

%% Before moving for the calculation of unmeasured variables we will assess
% if any of the measured variable is non-redundant, and therefore, cannot 
% be reconciled directly
% NOTE that this section of redundancy check has to run only once at the
% beginning, proper arrangement is made to ensure the same
if (redundancy_check_status == 0)
    for i = 1:length(right_second_part)
        if (right_second_part(i) == 0)
            nr_counter                = nr_counter + 1;         % counting the number of non-redundant measured variables starting from 0
            nr_positions(nr_counter)  = i;                      % recording the position of non-redundant measured variable in a vector
        end
    end
    if (nr_counter == 0)
        fprintf('\nNON-REDUNDANCY: All the measured variables are REDUNDANT\n\n')                                                                       % Test failed
    else
        nr_sym = M_sym(nr_positions);                   % creates a vector of symbolic variable from M_sym for non-red variables
        fprintf('\n\nNON-REDUNDANCY: There are %d measured variables out of %d that are non-redundant.\n',nr_counter,length(right_second_part))
        disp(['NON-REDUNDANCY: The NON-REDUNDANT variables so found are listed below.'])
        disp(nr_sym)                                    % will display the symbolic vector containing names of non-red. variables                                                                    % Test passed
    end  
    redundancy_check_status = 1;                        % this says that the on time redundancy check is performed
    if (length(non_redundants_equations) ~= nr_counter) % non-redundant variables need to be calulcated from equality equations. Here we check if provided number of equations equal to the number of non-redundant measured variables
        SAE_M = 0; SAE_U = 0;                           % allocating trial values just to end the loop without error
        return
    end
end

%% Here we are reconciling the non-redundant measured variables based on redundant variables
if (nr_counter > 0)                                          % if there are non-redundant measured variables, then they need manual adjustment
    Uhat = -inv(Au'*Au)*Au'*(Am*Mhat-b);                     % getting calculated unmeasured variables from reconciled measured variables
    %Uhat       = -inv_R1* (Q1' * (Am * Mhat - b));          % getting calculated unmeasured variables from reconciled measured variables; this equation is valid when non-observable variables are present
    % Filling up symbolic variables with latest generated values for later
    % subsitutions in equations provided for non-redundant variables
    for i = 1:length_M
        evalc(sprintf('%s=Mhat(%d)',char(M_sym(i)),i));            
    end
    for i = 1:length_U
        evalc(sprintf('%s=Uhat(%d)',char(U_sym(i)),i));
    end
    Mhat(nr_positions) = double(subs(non_redundants_equations));           % replacing by the correct values obtained from additional equations 
end

%% reCalculation for unmeasured variables 
Uhat        = -inv(Au'*Au)*Au'*(Am*Mhat-b);              % getting calculated unmeasured variables from reconciled measured variables
%Uhat       = -inv_R1* (Q1' * (Am * Mhat - b));          % getting calculated unmeasured variables from reconciled measured variables; this equation is valid when non-observable variables are present

%% Some calculations for convergence check
Obj_fun_current = (Mhat' - M_0)*W*(Mhat' - M_0)';  % current Obj_fun value
SAE_M       = sum(abs(Mhat'-M));                   % difference of last consequetive y values for convergence check
SAE_U       = sum(abs(Uhat'-U));                   % difference of last consequetive z values for convergence check
mean_abs_f  = mean(abs(f_MU));                     % average value of equality contraint function or fun_MU or f_MU
std_f_MU    = std(abs(f_MU));                      % standard deviation of f_MU value; a value of zero indicate proper fitting of iterative results into equality constraints

%% Filling the result matrix with obtained values
result(n_rows_result+1,1:length_M)                      = Mhat';
result(n_rows_result+1,length_M+1:length_M+length_U)    = Uhat';
result(n_rows_result+1,length_M+length_U+1)             = Obj_fun_current;                                      
result(n_rows_result+1,length_M+length_U+1+1)           = SAE_M;            % difference of last consequetive y values for convergence check
result(n_rows_result+1,length_M+length_U+1+1+1)         = SAE_U;            % difference of last consequetive z values for convergence check
result(n_rows_result+1,length_M+length_U+1+1+1+1)       = mean_abs_f;       % average value of equality contraint function or fun_MU or f_MU
result(n_rows_result+1,length_M+length_U+1+1+1+1+1)     = std_f_MU;         % standard deviation of f_MU value; a value of zero indicate proper fitting of iterative results into equality constraints

%% Plotting
iteration_vector = (2:(n_rows_result+1))-1;

subplot(3,2,[1 2])
bar([Mhat' Uhat'])
set(gca,'xtick',1:length([Mhat' Uhat']))
set(gca,'xticklabel',x_axis_marks_measured_unmeasured,'fontsize',labelfontsize)
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',axisfontsize)
% y_max = 1.1*max(M_0);
% ylim([-y_max y_max])   
grid on
grid MINOR
title('Instantaneous value of variables','fontsize',titlefontsize)
%ylabel('values','fontsize',labelfontsize)

subplot(3,2,3)
plot(iteration_vector,result(2:(n_rows_result+1),length_M+length_U+1),'--gs',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0.5,0.5,0.5])
set(gca,'xtick',1:iteration_vector(end))
set(gca,'XTickLabelRotation',45)
value_vector = result(2:(n_rows_result+1),length_M+length_U+1);
lgnd_Obj_fun = legend(sprintf('CurrentValue = %0.3e',value_vector(end)));
set(lgnd_Obj_fun,'FontSize',14);
title('Objective function value','fontsize',titlefontsize)
grid on
grid MINOR
set(gca,'FontSize',axisfontsize)
xlabel('iteration','fontsize',labelfontsize)
%ylabel('Obj_fun','fontsize',labelfontsize)

subplot(3,2,4)
plot([1:length(f_MU)],f_MU,'--gs',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0.5,0.5,0.5])
set(gca,'xtick',1:length(f_MU))
set(gca,'XTickLabelRotation',45)
title('Instantaneous equality function values (f)','fontsize',titlefontsize)
set(gca,'FontSize',axisfontsize)
grid on
grid MINOR
xlabel('Equality functions in sequence','fontsize',labelfontsize)
%ylabel('$f$', 'interpreter','latex','fontsize',labelfontsize)

subplot(3,2,5)
hold on
plot(iteration_vector,result(2:(n_rows_result+1),length_M+length_U+1+1+1+1+1),'--gs',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0.5,0.5,0.5])
plot(iteration_vector,result(2:(n_rows_result+1),length_M+length_U+1+1+1+1),'--rs',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[0.5,0.5,0.5])
set(gca,'xtick',1:iteration_vector(end))
set(gca,'XTickLabelRotation',45)
set(gca, 'YScale', 'log')
value_vector_sd_f   = result(2:(n_rows_result+1),length_M+length_U+1+1+1+1+1);
value_vector_mean_f = result(2:(n_rows_result+1),length_M+length_U+1+1+1+1);
lgnd_sd_mean_f = legend(sprintf('CurrentValue (\\sigma) = %0.3e',value_vector_sd_f(end)), sprintf('CurrentValue (\\mu) = %0.3e',value_vector_mean_f(end)));
set(lgnd_sd_mean_f,'FontSize',14);
title({'Standard deviation and mean of'; 'absolute values of equality functions'},'fontsize',titlefontsize)
set(gca,'FontSize',axisfontsize)
xlabel('iteration','fontsize',labelfontsize)
grid on
grid MINOR
%ylabel('\sigma','fontsize',labelfontsize)

subplot(3,2,6)
hold on 
plot(iteration_vector,result(2:(n_rows_result+1),length_M+length_U+1+1),'--gs',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0.5,0.5,0.5])
plot(iteration_vector,result(2:(n_rows_result+1),length_M+length_U+1+1+1),'--rs',...
        'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[0.5,0.5,0.5])
set(gca,'xtick',1:iteration_vector(end))
set(gca,'XTickLabelRotation',45)
set(gca, 'YScale', 'log')
value_vector_measured   = result(2:(n_rows_result+1),length_M+length_U+1+1);
value_vector_unmeasured = result(2:(n_rows_result+1),length_M+length_U+1+1+1);
lgnd_SAE = legend(sprintf('CurrentValue (SAE(M)) = %0.3e',value_vector_measured(end)), sprintf('CurrentValue (SAE(U)) = %0.3e',value_vector_unmeasured(end)));
set(lgnd_SAE,'FontSize',14);
title({'Convergence plot for'; 'measured and unmeasured variables'},'fontsize',titlefontsize)
set(gca,'FontSize',axisfontsize)
xlabel('iteration','fontsize',labelfontsize)
grid on
grid MINOR 
%ylabel('residual','fontsize',labelfontsize)

drawnow

%% Display imp values
if(iteration_vector(end) == 1)                      % shall print title only once
    %output to screen
    fprintf('Iterations Convergence Progress \n')               
    fprintf('iteration	mean(|f|)     Standard_deviation(f) SAE_M       SAE_U       Obj_fun \n')                        % iteration(end) returns the last values of iteration_vector 

    % deleting the existing file
    if exist('outputData.csv', 'file')==2
        delete('outputData.csv');
    end
    %output to file
    title_bar = sprintf('iteration    mean(|f|)     std(f)      SAE_M           SAE_U           Obj_fun         %s         ',strjoin(x_axis_marks_measured_unmeasured));
    dlmwrite('outputData.csv',title_bar,'-append','delimiter','')
end
%output to screen
fprintf( 'iter:   %d       %10f   %10f          %10f      %10f        %10f \n',iteration_vector(end), mean_abs_f, std_f_MU, SAE_M, SAE_U, Obj_fun_current)                       % iteration(end) returns the last values of iteration_vector

%output to file
outputData = sprintf('%.10f ', [iteration_vector(end) mean_abs_f std_f_MU,SAE_M SAE_U Obj_fun_current Mhat' Uhat']);
dlmwrite('outputData.csv',outputData,'-append','delimiter','')            % dlmwrite was used because fprintf write does not work in appending mode properly