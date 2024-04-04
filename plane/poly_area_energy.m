function [energy,A] = poly_area_energy( x, y) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the centroid coordinates of a Voronoi region as
% well as its individual Voronoi energy
%
% INPUT: - x,y are the coordinates of the vertices of the Voronoi cell
%        - genx, geny are the coordinates of the generator of the cell.
%    
% OUTPUT: - centroid is a list containg the coordinates of centroid
%           of the cell (x_cen and y_cen).
%         - energy is the individual energy of the Voronoi cell
%
% Initially written by H.J. Sommer III - 16.12.09 - tested under MATLAB v9.0
% Modified by Ivan Gonzalez @McGill University 
% Last update: May 17th of 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if inputs are same size
if ~isequal( size(x), size(y) )
  error( 'X and Y must be the same size');
end
 
% We shift the vertices by [genx;geny] to compute quantities with respect
% to the generator and not the origin of the domain.
x = x - mean(x);
y = y - mean(y);
  
% summations for CCW boundary. Direct application of Green's theorem to
% compute quntities over the polygon
xp = x( [2:end 1] );
yp = y( [2:end 1] );
a = x.*yp - xp.*y;

b=sqrt((x-xp).^2+(y-yp).^2);
P=sum(b); % positive perimeter
 
A = sum( a ) /2;
xc = sum( (x+xp).*a  ) /6/A;
yc = sum( (y+yp).*a  ) /6/A;
Ixx = sum( (y.*y +y.*yp + yp.*yp).*a  ) /12;
Iyy = sum( (x.*x +x.*xp + xp.*xp).*a  ) /12;

A=abs(A);

% In case we have a CW order we change the sigh

Ixx=abs(Ixx); Iyy=abs(Iyy);
energy=Ixx+Iyy;
 

end
