function [c,ceq] = DR_SQP_Constraints_VandZn(x)

c = [];

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

% M_Qs    =   [Qs1        Qs14        Qs16        Qs30        Qs33];
% x       =   [x(1)       x(2)        x(3)        x(4)        x(5)];
%U_Qs     =   [Qs7        Qs12        Qs13        Qs15        Qs18        Qs26        Qs31        Qs36        Qs37];
% x       =   [x(6)       x(7)        x(8)        x(9)        x(10)       x(11)       x(12)       x(13)       x(14)];
% M_Ql    =   [Ql3        Ql4         Ql8         Ql9         Ql11        Ql17        Ql20        Ql28        Ql29        Ql32        Ql35];
% x       =   [x(15)      x(16)       x(17)       x(18)       x(19)       x(20)       x(21)       x(22)       x(23)       x(24)       x(25)];
% U_Ql    =   [Ql19       Ql21        Ql27        Ql34        Ql38];
% x       =   [x(26)      x(27)       x(28)       x(29)       x(30)];
% M_Cw    =   [Cw7        Cw12        Cw13        Cw14        Cw15        Cw16        Cw18        Cw26        Cw31        Cw33];
% x       =   [x(31)      x(32)       x(33)       x(34)       x(35)       x(36)       x(37)       x(38)       x(39)       x(40)];
% M_Cs    =   [CS1        CS14        CS16        CS30        CS33        CS7         CS12        CS13        CS15        CS18        CS26        CS31        CS36        CS37];
% x       =   [x(41)      x(42)       x(43)       x(44)       x(45)       x(46)       x(47)       x(48)       x(49)       x(50)       x(51)       x(52)       x(53)       x(54)];
% M_Cl    =   [CL4        CL7         CL8         CL9         CL11        CL12        CL13        CL14        CL15        CL16        CL18        CL19        CL20        CL21        CL26        CL27      CL29        CL31        CL33        CL34];
% x       =   [x(55)      x(56)       x(57)       x(58)       x(59)       x(60)       x(61)       x(62)       x(63)       x(64)       x(65)       x(66)       x(67)       x(68)       x(69)       x(70)     x(71)       x(72)       x(73)       x(74)];

%% Commmon equality constraints
Common_equalities = [x(1)-x(6);  % solid flow rate balancing
    x(7)+x(2)-x(9);
    x(9)+x(3)-x(10);
    x(4)+x(10)-x(2)-x(3)-x(12);
    x(8)+x(12)-x(5);
    x(11)-x(13);
    x(13)-x(14);

    % Solution mass balance (Balancing the liquid only flow rates
    x(1)*(1/x(115)-1)+ x(15)+x(16)-x(6)*(1/x(31)-1)+x(17);
    x(7)*(1/x(32)-1)+x(2)*(1/x(34)-1)-x(9)*(1/x(35)-1);
    x(9)*(1/x(35)-1)+x(3)*(1/x(36)-1)+x(20)-x(10)*(1/x(37)-1)-x(26);
    x(10)*(1/x(37)-1)- x(2)*(1/x(34)-1)-x(3)*(1/x(36)-1)-x(12)*(1/x(39)-1)+x(29)+x(117);
    x(8)*(1/x(33)-1)- (x(19)+x(117))+x(12)*(1/x(39)-1)+x(24)-x(5)*(1/x(40)-1);   
    x(26)+x(21)-x(27);
    x(23)-x(18)-x(16)-x(30);  
    x(28)-x(17)+x(22)-x(23);
    x(11)*(1/x(38)-1)-x(29)+x(25)-x(13)*(1/x(116)-1);

    % Total balance at digestion-sand separation and precipitation circuit
    % Normal
    % x(6)/x(31)+x(18)+x(19)-x(7)/x(32)-x(8)/x(33);
    % DigestionSandSeperated
    x(6)+x(6)*((1/x(31))-1)+x(18)-x(75)-x(75)*((1/x(76))-1);
    x(75)+x(75)*((1/x(76))-1)+x(19)-x(7)-x(7)*((1/x(32))-1)-x(8)-x(8)*((1/x(33))-1);
    x(27)-x(11)-x(11)*(1/x(38)-1)-x(28);
    x(34)-x(36)];

%% Vanadium species balance
Zn_equalities = [x(1)*x(41)+x(16)*x(55)-x(6)*x(46)-x(6)*(1/x(31)-1)*x(56)+x(17)*x(57);
    
    % Total balance at digestion-sand separation and precipitation circuit
    % Normal or 
    % x(6)*x(46)-x(6)*(1/x(31)-1)*x(56)+x(18)*x(58)+x(19)*x(59)-x(7)*x(47)-x(7)*(1/x(32)-1)*x(60)-x(8)*x(48)-x(8)*(1/x(33)-1)*x(61);
    % DigestionSandSeperated
    x(6)*x(46)+x(6)*((1/x(31))-1)*x(56)+x(18)*x(58)-x(75)*x(77)-x(75)*((1/x(76))-1)*x(78);
    x(75)*x(77)+x(75)*((1/x(76))-1)*x(78)-x(7)*x(47)-x(7)*((1/x(32)-1))*x(60)-x(8)*x(48)-x(8)*((1/x(33)-1))*x(61)+ x(19)*x(59);
    x(7)*x(47)+x(7)*(1/x(32)-1)*x(60)+x(2)*x(42)+x(2)*(1/x(34)-1)*x(62)-x(9)*x(49)-x(9)*(1/x(35)-1)*x(63);    
    x(9)*x(49)+x(9)*(1/x(35)-1)*x(63)+x(3)*x(43)+x(3)*(1/x(36)-1)*x(64)-x(10)*x(50)-x(10)*(1/x(37)-1)*x(65)-x(26)*x(66);    
    x(26)*x(66)+x(21)*x(67)-x(27)*x(68);   
    x(28)*x(70)-x(17)*x(57)-x(23)*x(71);    
    x(23)*x(71)-x(16)*x(55)-x(18)*x(58);   
    x(10)*x(50)+x(10)*(1/x(37)-1)*x(65)-x(2)*x(42)-x(2)*(1/x(34)-1)*x(62)-x(3)*x(43)-x(3)*(1/x(36)-1)*x(64)-x(12)*x(52)-x(12)*(1/x(39)-1)*x(72)+x(29)*x(74)+x(117)*x(44);    
    x(8)*x(48)+x(8)*(1/x(33)-1)*x(61)+x(12)*x(52)+x(12)*(1/x(39)-1)*x(72)-x(19)*x(59)-x(117)*x(44)-x(5)*x(45)-x(5)*(1/x(40)-1)*x(73);    
    x(11)*x(51)+x(11)*(1/x(38)-1)*x(69)-x(29)*x(74)-x(13)*x(53);    
    x(13)*x(53)-x(14)*x(54);
    x(27)*x(68)-x(11)*x(51)-x(11)*(1/x(38)-1)*x(69)-x(28)*x(70);
    
    % Additional supporting constraints
    x(55)-x(58);
    x(57)-x(71);
    x(42)-x(43);
    x(62)-x(64);
    x(44)-x(59);
    x(53)-0.6537*x(54)]; 

%% Zinc species balance
V_equalities = [x(1)*x(79)+x(16)*x(93)-x(6)*x(84)-x(6)*(1/x(31)-1)*x(94)+x(17)*x(95);
    
    % Total balance at digestion-sand separation and precipitation circuit
    % Normal or 
    % x(6)*x(46)-x(6)*(1/x(31)-1)*x(56)+x(18)*x(58)+x(19)*x(59)-x(7)*x(47)-x(7)*(1/x(32)-1)*x(60)-x(8)*x(48)-x(8)*(1/x(33)-1)*x(61);
    % DigestionSandSeperated
    x(6)*x(84)+x(6)*((1/x(31))-1)*x(94)+x(18)*x(96)-x(75)*x(113)-x(75)*((1/x(76))-1)*x(114);
    x(75)*x(113)+x(75)*((1/x(76))-1)*x(114)-x(7)*x(85)-x(7)*((1/x(32)-1))*x(98)-x(8)*x(86)-x(8)*((1/x(33)-1))*x(99)+x(19)*x(97);
    x(7)*x(85)+x(7)*(1/x(32)-1)*x(98)+x(2)*x(80)+x(2)*(1/x(34)-1)*x(100)-x(9)*x(87)-x(9)*(1/x(35)-1)*x(101);    
    x(9)*x(87)+x(9)*(1/x(35)-1)*x(101)+x(3)*x(81)+x(3)*(1/x(36)-1)*x(102)-x(10)*x(88)-x(10)*(1/x(37)-1)*x(103)-x(26)*x(104);    
    x(26)*x(104)+x(21)*x(105)-x(27)*x(106);   
    x(28)*x(108)-x(17)*x(95)-x(23)*x(109);    
    x(23)*x(109)-x(16)*x(93)-x(18)*x(96);   
    x(10)*x(88)+x(10)*(1/(37)-1)*x(103)-x(2)*x(80)-x(2)*(1/x(34)-1)*x(100)-x(3)*x(81)-x(3)*(1/x(36)-1)*x(102)-x(12)*x(90)-x(12)*(1/x(39)-1)*x(110)+x(29)*x(112)+x(117)*x(82);    
    x(8)*x(86)+x(8)*(1/x(33)-1)*x(99)+x(12)*x(90)+x(12)*(1/x(39)-1)*x(110)-x(19)*x(97)-x(117)*x(44)-x(5)*x(83)-x(5)*(1/x(40)-1)*x(111);    
    x(11)*x(89)+x(11)*(1/x(38)-1)*x(107)-x(29)*x(112)-x(13)*x(91);    
    x(13)*x(91)-x(14)*x(92);
    x(27)*x(106)-x(11)*x(89)-x(11)*(1/x(38)-1)*x(107)-x(28)*x(108);
    
    % Additional supporting constraints
    x(93)-x(96);
    x(95)-x(109);
    x(80)-x(81);
    x(100)-x(102);
    x(97)-x(82);
    x(91)-0.6537*x(92)]; 

%% All equality constraints
ceq = [Common_equalities;
        Zn_equalities;
        V_equalities];

%% Flow rate equalities
% Common_equalities = [Qs1-Qs7;  % solid flow rate balancing
%     Qs12+Qs14-Qs15;
%     Qs15+Qs16-Qs18;
%     Qs14+Qs16-Qs18-Qs30+Qs31;
%     Qs13+Qs31-Qs33;
%     Qs26-Qs36;
%     Qs36-Qs37;
%     Qs1+Qs30-Qs26-Qs33;
%
%     % Solution mass balance (Balancing the liquid only flow rates)
%     Ql3+Ql4-Qs7*(1/Cw7-1)+Ql8;
%     Qs12*(1/Cw12-1)+Qs14*(1/Cw14-1)-Qs15*(1/Cw15-1);
%     Qs15*(1/Cw15-1)+Qs16*(1/Cw16-1)+Ql17-Qs18*(1/Cw18-1)-Ql19;
%     Qs14*(1/Cw14-1)+Qs16*(1/Cw16-1)-Qs18*(1/Cw18-1)+Qs31*(1/Cw31-1)-Ql34;
%     Ql11-Qs13*(1/Cw13-1)-Qs31*(1/Cw31-1)-Ql32+Qs33*(1/Cw33-1);
%     Ql19+Ql20-Ql21;
%     Ql4+Ql9-Ql29+Ql38;
%     (Ql27-Ql8)+Ql28-Ql29;
%     Qs26*(1/Cw26-1)-Ql34+Ql35;
%  
%     % Total balance for digestion-sand separation and precipitation
%     circuit
%     Qs7/Cw7+Ql9+Ql11-Qs12/Cw12-Qs13/Cw13;
%     Ql21-Qs26-Qs26*(1/Cw26-1)-Ql27;
%
%     % Overall solid-liquid balance
%     Ql3+Ql28+Ql17+Ql20+Ql32+Ql35-Qs33*(1/Cw33-1)-Ql38;
%
%     % Additional supporting constraints
%     Cw14-Cw16];

%% Vanadium species balance
% V_equalities = [Qs1*CsV1+Ql4*ClV4-Qs7*CsV7-Qs7*(1/Cw7-1)*ClV7+Ql8*ClV8;    
%     Qs7*CsV7+Qs7*(1/Cw7-1)*ClV7+Ql9*ClV9+Ql11*ClV11-Qs12*CsV12-Qs12*(1/Cw12-1)*ClV12-Qs13*CsV13-Qs13*(1/Cw13-1)*ClV13;    
%     Qs12*CsV12+Qs12*(1/Cw12-1)*ClV12+Qs14*CsV14+Qs14*(1/Cw14-1)*ClV14-Qs15*CsV15-Qs15*(1/Cw15-1)*ClV15;
%     Qs15*CsV15+Qs15*(1/Cw15-1)*ClV15+Qs16*CsV16+Qs16*(1/Cw16-1)*ClV16-Qs18*CsV18-Qs18*(1/Cw18-1)*ClV18-Ql19*ClV19;
%     Ql19*ClV19+Ql20*ClV20-Ql21*ClV21;
%     (Ql27-Ql8)*ClV8-Ql29*ClV29;
%     Ql4*ClV4+Ql9*ClV9-Ql29*ClV29;
%     Qs14*CsV14+Qs14*(1/Cw14-1)*ClV14+Qs16*CsV16+Qs16*(1/Cw16-1)*ClV16-Qs18*CsV18-Qs18*(1/Cw18-1)*ClV18+Qs31*CsV31+Qs31*(1/Cw31-1)*ClV31-Ql34*ClV34+Qs30*CsV30;
%     Ql11*ClV11-Qs13*CsV13-Qs13*(1/Cw13-1)*ClV13-Qs31*CsV31-Qs31*(1/Cw31-1)*ClV31+Qs33*CsV33+Qs33*(1/Cw33-1)*ClV33;
%     Qs26*CsV26+Qs26*(1/Cw26-1)*ClV26-Ql34*ClV34-Qs36*CsV36;
%     Qs36*CsV36-Qs37*CsV37;
%     Ql21*ClV21-Qs26*CsV26-Qs26*(1/Cw26-1)*ClV26-Ql27*ClV27;
%     % Additional supporting constraints
%     ClV4-ClV9;
%     ClV8-ClV27;
%     CsV14-CsV16;
%     ClV14-ClV16; 
%     CsV36-CsV37]; 
                    
%% Valuable Zinc content balance
% Zn_equalities = [Qs1*CsZn1+Ql4*ClZn4-Qs7*CsZn7-Qs7*(1/Cw7-1)*ClZn7+Ql8*ClZn8;
%     Qs7*CsZn7+Qs7*(1/Cw7-1)*ClZn7+Ql9*ClZn9+Ql11*ClZn11-Qs12*CsZn12-Qs12*(1/Cw12-1)*ClZn12-Qs13*CsZn13-Qs13*(1/Cw13-1)*ClZn13;    
%     Qs12*CsZn12+Qs12*(1/Cw12-1)*ClZn12+Qs14*CsZn14+Qs14*(1/Cw14-1)*ClZn14-Qs15*CsZn15-Qs15*(1/Cw15-1)*ClZn15;
%     Qs15*CsZn15+Qs15*(1/Cw15-1)*ClZn15+Qs16*CsZn16+Qs16*(1/Cw16-1)*ClZn16-Qs18*CsZn18-Qs18*(1/Cw18-1)*ClZn18-Ql19*ClZn19;
%     Ql19*ClZn19+Ql20*ClZn20-Ql21*ClZn21;
%     (Ql27-Ql8)*ClZn8-Ql29*ClZn29;
%     Ql4*ClZn4+Ql9*ClZn9-Ql29*ClZn29;
%     Qs14*CsZn14+Qs14*(1/Cw14-1)*ClZn14+Qs16*CsZn16+Qs16*(1/Cw16-1)*ClZn16-Qs18*CsZn18-Qs18*(1/Cw18-1)*ClZn18+Qs31*CsZn31+Qs31*(1/Cw31-1)*ClZn31-Ql34*ClZn34+Qs30*CsZn30;
%     Ql11*ClZn11-Qs13*CsZn13-Qs13*(1/Cw13-1)*ClZn13-Qs31*CsZn31-Qs31*(1/Cw31-1)*ClZn31+Qs33*CsZn33+Qs33*(1/Cw33-1)*ClZn33;
%     Qs26*CsZn26+Qs26*(1/Cw26-1)*ClZn26-Ql34*ClZn34-Qs36*CsZn36;
%     Qs36*CsZn36-Qs37*CsZn37;
%     Ql21*ClZn21-Qs26*CsZn26-Qs26*(1/Cw26-1)*ClZn26-Ql27*ClZn27;
%     
%     % Additional supporting constraints
%     ClZn4-ClZn9;
%     ClZn8-ClZn27;
%     CsZn14-CsZn16;
%     ClZn14-ClZn16;
%     CsZn36-CsZn37];
