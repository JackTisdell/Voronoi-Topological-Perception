function value = arc_cosine(cc)

%*******************************************************************************
%
%% ARC_COSINE computes the arc cosine function, with argument truncation--> to avoid complex numbers.

%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%    Author: John Burkardt
%    Modified: 28 January 2005
%
%  Parameters:
%
%    Input, real C, the argument.
%    Output, real VALUE, an angle whose cosine is C.
%
  cc2 = cc;
  cc2 = max ( cc2, -1.0 );
  cc2 = min ( cc2, +1.0 );

  value = acos ( cc2 );

  return
end