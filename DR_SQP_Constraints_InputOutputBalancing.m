function [c,ceq] = DR_SQP_Constraints_GrossErrorDemo(x)

c = [];

%% Commmon equality constraints
% Common_equalities = [x(1)-x(2)-x(3);]; % Collection 5
 Common_equalities = [x(1)+x(2)+x(3)-x(4)+x(5)+x(6)+x(7)+x(8)+x(9)-x(10)-x(10)/x(12)-x(11);]; % Collection 6

%% Vanadium species balance
Zn_equalities = []; 

%% Zinc species balance
V_equalities = []; 

%% All equality constraints
ceq = [Common_equalities;
        Zn_equalities;
        V_equalities];