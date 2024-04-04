function A = poly_area(x,y) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the centroid coordinates of a Voronoi region as
% well as its individual Voronoi energy
%
% INPUT: - x,y are the coordinates of the vertices of the polygon
%        
%    
% OUTPUT: - A is area of the polygon
%         
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
 
A = sum( a ) /2;
A=abs(A); 

end
