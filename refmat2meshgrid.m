% Take a 3x2 Matlab reference matrix (R) and finds the meshgrid x,y
% coordinates
% 
% RELEASE NOTES
%   Written by Mark Raleigh (mraleig1@uw.edu), January 2012
%   Revised by Mark Raleigh (April 2014) to address potential bug
% 
% SYNTAX
%   [X,Y]= refmat2meshgrid(R,rows,cols)
% 
% INPUTS
%   R = 3x2 Matlab reference matrix (see makerefmat)
%   rows = number of rows in the output
%   cols = number of columns in the output
% 
% OUTPUTS
%   X = x-coordinate matrix
%   Y = y-coordinate matrix

function [X,Y]= refmat2meshgrid(R,rows,cols)

dx = R(2,1);
dy = R(1,2);

ULcenter_X = R(3,1)+dx;
URcenter_X = R(3,1)+(dx*cols);
% x = ULcenter_X:dx:URcenter_X;
x = linspace(ULcenter_X, URcenter_X, cols);

ULcenter_Y = R(3,2)+dy;
LLcenter_Y = R(3,2)+(dy*rows);
% y = ULcenter_Y:dy:LLcenter_Y;
y = linspace(ULcenter_Y, LLcenter_Y, rows);


[X,Y] = meshgrid(x,y);