# DR
Data Reconciliation MATLAB program by Successive Linearization (SL) and Sequential Quadratic Programming (SQP)

DR_SQP_Main 
(filename,SQPSheetName,x0range,SDrange,~,lbrange,ubrange,equality_write_range,x_writerange)
●	The primary code file for the Data Reconciliation

DR_SQP_Constraints_VorZn (x)
●	Contains the mass balance equations for either Vanadium or Zinc component balance, either of them can be utilized

DR_SQP_Constraints_VandZn (x)
●	Contains the mass balance equations for simultaneous Vanadium and Zinc component balance, both of them to be utilized

DR_SQP_Constraints_InputOutputBalancing(x)
●	Contains mass balance equations for only flow and no component balance 

Function application (app,event)
●	The whole application is defined inside this function which further used all above defined functions to calculate the reconciled data.
●	This is specifically the user interface function.(front-end)
