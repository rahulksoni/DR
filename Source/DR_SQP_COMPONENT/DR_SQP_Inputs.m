%% Inputs for SQP program for Data Reconciliation
% xlfilename          = '../../../NALCOProjectExcelDoNotRelocate';
% SQPSheetName        = 'SQPInputIndividualCollection';
% x0range             = 'p7:p80';
% SDrange             = 'f7:f80';
% lbrange             = 'c7:c80';
% ubrange             = 'd7:d80';
% x_writerange        = 'r7:r80';
% equality_write_range= 'an7:an44';

% xlfilename          = '../../../DRExcelFiles/ReconcilingZinc/V3_ZincResults';
xlfilename          = '../../../DRExcelFiles/Reconciling_VZnCombined/V3_V_Zn_Combined_ChangesAfterNALCOSuggestion';
SQPSheetName        = 'SQPInputIndividualCollection';
x_writerange        = 'r11:r127';
equality_write_range= 'an11:an127';

collection = 2;
% Collection 2, 3, 4 are actual collection of NALCO plant
% Collection 5 and above correspond to Input Output balancing sections

if (collection == 2)   
    SDrange             = 'e11:e127';
    lbrange             = 'k11:k127';
    ubrange             = 'l11:l127';
    x0range             = 'm11:m127';
elseif (collection == 3)  
    SDrange             = 'e11:e127';
    lbrange             = 'z11:z127';
    ubrange             = 'aa11:aa127';
    x0range             = 'ab11:ab127';
elseif (collection == 4)
    SDrange             = 'e11:e127';
    lbrange             = 'ap11:ap127';
    ubrange             = 'aq11:aq127';
    x0range             = 'ar11:ar127';
elseif (collection == 5)  
    SDrange             = 'e142:e144';
    lbrange             = 'k142:k144';
    ubrange             = 'l142:l144';
    x0range             = 'm142:m144';
elseif (collection == 6)  
    SDrange             = 'e147:e158';
    lbrange             = 'k147:k158';
    ubrange             = 'l147:l158';
    x0range             = 'm147:m158';
end