function [x0,final_values,x,ceq_x] = DR_SQP_Main(xlfilename,SQPSheetName,x0range,SDrange,~,lbrange,ubrange,equality_write_range,x_writerange)
x0=[];
x=[];
final_values=[];
ceq_x=[];
% % FileName = inputs_Attempt2_NALCOBothCase.m
%clear all                   % clear the variables data
close all                   % closes the opened figure windows
format SHORTENG
    
%% Input data/materi f file name for input data updated in file inputs.m is imported in program through command eval
%input_file_name = 'DR_SQP_Inputs';
%input_file_name = 'DR_SQP_Inputs_GoldCase';
%eval(input_file_name);  

%% Data extraction part
x0          = xlsread(xlfilename,SQPSheetName,x0range)';   % extracting initial values from excel file
SD          = xlsread(xlfilename,SQPSheetName,SDrange)';   % extracting standard deviation values from excel file

%% Preparing co-variance matrix
%V           = power(SD,2);                                      % Variance vector for the measured variables
V           = SD.^2;
W           = diag(1./V);                                       % Co-variance matrix of inverse variance values 

%% Data checkup
zero_SD_value_flag = 0;                                         % Shall turn to 1 if any of the SD value is zero, SD values cannot be zero as there inverse has to be taken in the W = diag(1./V);
length_SD   = length(SD);
for i=1:length_SD
    if (SD(i) == 0) 
       zero_SD_value_flag = 1;
       fprintf('The SD value at position %d is 0. Consider converting it to a small value',i)
       return
    end
end

%% The run can be done for either V or Zn concerntraion or for the combination of both of them
% Lets user tell which one s/he wants
fprintf('Whether V or Zn or Both of them are required to be run? \n');
prompt = 'Enter 0 for V or Zn, otherwise 1 for V-Zn combined run or 2 for Input Output Balancing? ';singleOrBoth = 1;
     
%% Optimization part
%% Objective Function
% fun       = @(x)((2*(x-x0)*W)*(x-x0)'+(x-x0)*2*W*(x-x0)');    % the primary objective function
fun         = @(x)((x-x0)*W*(x-x0)');                           % the primary objective function
A           = [];                                               % coefficients to be written in A*x<b
b           = [];                                               % inequality bound to be written in A*x<b
Aeq         = [];                                               % coefficient to be written in Aeq*x=beq
beq         = [];                                               % equality value to be written in Aeq*x=beq
lb          = xlsread(xlfilename,SQPSheetName,lbrange)';        % lower bound on the variables
ub          = xlsread(xlfilename,SQPSheetName,ubrange)';        % upper bound on the variables
if (singleOrBoth == 1)
    nonlcon     = @DR_SQP_Constraints_VorZn;                    % a separate files to be called containing non-linear equality and inequality constraints
elseif (singleOrBoth == 1)
    nonlcon     = @DR_SQP_Constraints_VandZn;                   % a separate files to be called containing non-linear equality and inequality constraints
elseif (singleOrBoth == 2)
    nonlcon     = @DR_SQP_Constraints_InputOutputBalancing;           % a separate files to be called containing non-linear equality and inequality constraints
end
% nonlcon     = @DR_SQP_Constraints_GoldCase;                              % a separate files to be called containing non-linear equality and inequality constraints
options     = optimset('Display','iter','FunValCheck','on','Algorithm','sqp', 'MaxIter', 10000, 'MaxFunEvals', 100000,'TolFun',1e-2,'TolCon',1e-2,'TolX',1e-10);    % options for optimization with sequential quadratic (sqp) algorithm 
% f(x) is the value of objective function at current iteration.
% fval is the value of objective function at solution x.
% Feasibility (Ideal value = 0): Maximum constraint violation, where satisfied inequality constraints count as 0.
% Steplength: Multiplicative factor that scales the search direction.
% Norm of step: Size of the current step.
% First-order optimality: First-order optimality is a measure of how close a point x is to optimal. Most Optimization Toolbox? solvers use this measure, though it has different definitions for different algorithms. First-order optimality is a necessary condition, but it is not a sufficient condition. In other words:
%   The first-order optimality measure must be zero at a minimum.
%   A point with first-order optimality equal to zero is not necessarily a minimum. (More information at http://www.mathworks.com/help/optim/ug/first-order-optimality-measure.html)
% exitflag is the indicator that says why optimization has stopped. Visit https://www.mathworks.com/help/optim/ug/fmincon.html to see the meaning of different exitflags. 
% OptimalityTolerance (TolFun): Termination tolerance on the first-order optimality, a positive scalar. The default is 1e-6. See First-Order Optimality Measure. For optimset, the name is TolFun. See Current and Legacy Option Name Tables. The OptimalityTolerance tolerance relates to the first-order optimality measure. Typically, if the first-order optimality measure is less than OptimalityTolerance, solver iterations end.
% ConstraintTolerance (TolCon): Tolerance on the constraint violation, a positive scalar. The default is 1e-6. See Tolerances and Stopping Criteria. For optimset, the name is TolCon. See Current and Legacy Option Name Tables.
% StepTolerance	      (TolX)  : Termination tolerance on x, a positive scalar. The default value for all algorithms except 'interior-point' is 1e-6; for the 'interior-point' algorithm, the default is 1e-10. See Tolerances and Stopping Criteria. For optimset, the name is TolX. See Current and Legacy Option Name Tables.
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options); % final optimization command to predict the optimum values x

%% Calculation of objective function and Constraints at initial point and solution
fun_x0      = fun(x0);
fun_x       = fun(x);
if (singleOrBoth == 0)
    [c,ceq_x0]  = DR_SQP_Constraints_VorZn(x0);
    [c,ceq_x]   = DR_SQP_Constraints_VorZn(x);
elseif (singleOrBoth == 1)
    [c,ceq_x0]  = DR_SQP_Constraints_VandZn(x0);
    [~,ceq_x]   = DR_SQP_Constraints_VandZn(x); 
elseif (singleOrBoth == 2)
    [~,ceq_x0]  = DR_SQP_Constraints_InputOutputBalancing(x0);
    [~,ceq_x]   = DR_SQP_Constraints_InputOutputBalancing(x); 
end
% [c,ceq_x0]  = DR_SQP_Constraints_GoldCase(x0);
% [c,ceq_x]   = DR_SQP_Constraints_GoldCase(x);

fprintf('The initial and final value of variables with deviations are printed below.\n')
Variables_values_inital_final_deviations = [x0' x' (x'-x0')./x0'*100];

for i=1:length(ceq_x)
    if(abs(ceq_x(i,1))>0.1)
        ceq_not_satisfied(i,1) = 1;
    else
        ceq_not_satisfied(i,1) = 0;
    end
end
fprintf('The initial and final value of constraints functions along with non-satisfaction indicators are printed below.\n')
Constraint_function_inital_final_nonSatisfactionIndicator = [ceq_x0 ceq_x ceq_not_satisfied];

fprintf('Initial value of objective function = %.3f\n',fun_x0)
fprintf('Final value of objective function = %.3f\n',fun_x)

mean_of_absolute_of_equality_function_values = mean(abs(ceq_x))

%% Writing results into file with specified sheet name
x_write        = xlswrite(xlfilename,transpose(x),SQPSheetName,x_writerange)';
equality_write = xlswrite(xlfilename,ceq_x,SQPSheetName,equality_write_range)';
final_values = transpose(x)
equality_function_values = ceq_x
mean_absolute_equality_function_values = mean(abs(ceq_x))
max_absolute_equality_function_values = max(abs(ceq_x));

%% Setting up font sizes that are going to be used in plots
titlefontsize = 18;
axisfontsize  = 16;
labelfontsize = 18;

figure
subplot(2,2,[1,2])
bar(x)
title('Variables values at final iteration','fontsize',titlefontsize)
set(gca,'FontSize',axisfontsize)
xlabel('Variables','fontsize',labelfontsize)
grid on
grid MINOR 
ylabel('Varaiable values','fontsize',labelfontsize)

subplot(2,2,3)
bar((x-x0)./x0*100)
title('Deviation of variable values at final iteration','fontsize',titlefontsize)
set(gca,'FontSize',axisfontsize)
xlabel('Variables','fontsize',labelfontsize)
grid on
grid MINOR 
ylabel('Deviation %','fontsize',labelfontsize)

subplot(2,2,4)
bar(ceq_x)
title('Equality constraints values at final iteration','fontsize',titlefontsize)
set(gca,'FontSize',axisfontsize)
xlabel('Equatity constraints','fontsize',labelfontsize)
grid on
grid MINOR 
ylabel('','fontsize',labelfontsize)
