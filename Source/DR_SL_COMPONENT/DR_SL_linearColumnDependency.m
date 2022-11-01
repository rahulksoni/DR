function [n_ind_col,ind_col,A_reduced]=DR_SL_linearColumnDependency(A)

%% Based on https://math.stackexchange.com/questions/1179108/how-to-get-rid-of-linearly-dependent-columns-of-a-given-matrix-on-matlab
% Read about linear dependency in matrix at http://dankalman.net/AUhome/classes/classesS12/linalg/dependencies.pdf
% Theroy,
% the rank command roughly works in the following way: the matrix is converted to a reduced row echelon form, and then the number of rows that are not all equal 
% to zero are counted. Matlab uses a tolerance to determine what is equal to zero. If there is uncertainty in the numbers, you may have to define what zero is, 
% e.g. if the absolute value of a number is less than 1e-5, you may consider that close enough to be zero. The default tolerance is usually very small, of order 1e-15.
% Ref: http://matlab.cheme.cmu.edu/2011/08/02/determining-linear-independence-of-a-set-of-vectors/

%% Working principle
% There is a inbuild command in MATLAB 'rref(A)' which produces the reduced row echelon form of A.
% The result is the position of linearly independent columns. Note that it chooses the first column out of two dependent columns.
% Further, simple operation can result in reduced matrix with only independent columns.

%% Working with rows
% The matrix itself can be imported in transpose forms, and processed with
% removal of dependent columns and then back-transposed to original form

%% Codes
[~,ind_col] = rref(A);   
% Command rref(A) prepares the reduced echelon form of A with default tolerance criteria. 
% ~ symbol determines that we are not at all interested in row positions
% ind_col is a vector giving column positions for the independent columns
% ind_col stands independent-columns

A_reduced = A(:, ind_col); 
% Simply, selecting the columns of A that are independent and preparing a
% new reduced size matrix

n_ind_col = size(A_reduced,2);