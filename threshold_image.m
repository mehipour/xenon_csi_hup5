function [img_th,mask] = threshold_image(img,thresh_factor)
% syntax:
% [img_th,mask] = threshold_image(img,thresho)
% This function thresolds an image by threshold_factor times the maximum 
% value found in the img image. 
%
% 04/06/2014: Created by Mehrdad Pourfathi
%
% updated by MP on 04/15/2017
% ouput is also the mask.

% created masked pH map
th = thresh_factor*max(img(:));
img_th = zeros(size(img));
Ni = size(img);
mask = img_th;
for ii = 1:Ni(1);
    for jj = 1:Ni(2);
        if img(ii,jj) > th
            img_th(ii,jj) = img(ii,jj);
            mask(ii,jj) = 1;
        else
            img_th(ii,jj) = NaN;
        end
    end
end