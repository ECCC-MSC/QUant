function [Xi,Yi,Zi] = invdistgrid(X,Y,Z,dx,varargin)
%INVDISTGRID Simple, robust gridding using inverse-distance interpolation.
%   [Xi,Yi,Zi] = invdistgrid[X,Y,Z,dx] Grids values in vector Z with 
%       coordinates X and Y into an array with spacing dx using inverse-
%       distance weighting (default power =2). All data between grid points 
%       are used in weighting. If no data exists between points, a NaN is 
%       entered. X and Y must be vectors and Z must be of the same length. 
%       Z can be an m by n array, with each column a dataset, and Zi will 
%       be a three-dimnesional array with n layers.
%       
%   [Xi,Yi,Zi] = invdistgrid[X,Y,Z,dx,p] Replaces the default weighting
%       power 2 with the integer p. The higher the value of p, the more
%       weight will be given to closer data values.
% 
%   NOTE: You can use FILLNANS.m or INPAINT_NANS.m to fill in the NaN 
%   values in Zi.
%
%   Ian M. Howat, Applied Physics Lab, University of Washington
% 	ihowat@apl.washington.edu
%   Version 1: 09-Jul-2007 13:47:46


%weighting power
p = 2;
if nargin == 5;
    p = varargin{1};
end

%round the data locations to the interpolation grid spacing
XX = round(X./dx).*dx;
YY = round(Y./dx).*dx;

XYZ = sortrows([XX YY X Y Z],[2 1]);
ind = [0; XYZ(2:end,2) == XYZ(1:end-1,2) & XYZ(2:end,1) == XYZ(1:end-1,1); 0];
if sum(ind) > 0
    fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
    fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);

    for k = 1 : length(fs)
        D = repmat(sqrt((XYZ(fs(k),1) -  XYZ(fs(k):fe(k),3)).^2+...
            (XYZ(fs(k),2)-XYZ(fs(k):fe(k),4)).^2),[1,size(Z,2)]);
        D(D == 0) = .00001;
        XYZ(fe(k),5:end) = sum(XYZ(fs(k):fe(k),5:end)./(D.^p))./ sum(1./(D.^p));
    end
    XYZ = XYZ(~ind(2:end),:);
end

Xi = min(XYZ(:,1)):dx:max(XYZ(:,1));
Yi = fliplr(min(XYZ(:,2)):dx:max(XYZ(:,2)))';

c = (XYZ(:,1) -Xi(1)+dx)./dx;
r = (Yi(1)-XYZ(:,2)+dx)./dx;
I = index([r,c],length(Yi));

for k = 1:size(Z,2)
    tmp = ones(length(Yi),length(Xi)).*NaN;
    tmp(I) = XYZ(:,4+k);
    Zi(:,:,k) = tmp;
end

function I = index(ij,m)
% I = index([ij],[M]); Returns column-wise index equivalent I of position i j within 
% an M row array.

 i = ij(:,1);
 j = ij(:,2);
 
I = (j.* m) - m + i;

%return I as a column vector
I = reshape(I,[length(I) 1]);