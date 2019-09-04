function imgi = interpolate_image(img,int_fact)
% syntax:
% imgi = interpolate_image(img,int_fact)
% This function takes an image and interpolated it in both dimensions by a
% factor of int_fact.
%
% 04/04/2014: Created by Mehrdad Pourfathi
%
% Update (1) by Mehrdad Pourfathi
% int_fact is now an array of two. if it is only one element then
% automatically both parameters are set the same.
%
% Update (2) by Mehrdad Pourfathi on 6/16/17
% performs cubic interpolation instead of linear

if length(int_fact) < 2
    int_fact(2) = int_fact(1);
end;

Nx = size(img);
x1 = 1:Nx(1);
y1 = 1:Nx(2);
x2 = linspace(1,Nx(1),Nx(1)*int_fact(1));
y2 = linspace(1,Nx(2),Nx(2)*int_fact(2));
[X1,Y1] = meshgrid(x1,y1);
[X2,Y2] = meshgrid(x2,y2);
imgi1 = interp2(X1,Y1,img',X2,Y2,'cubic');
imgi = imgi1';

