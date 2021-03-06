function Chi2 = chistart (D,L,ahat,ncands,factor)
%CHISTART: Computes the initial size of the search ellipsoid
%
%    Chi2 = chistart (D,L,ahat,ncands,factor)
%
% This routine computes or approximates the initial size of the search
% ellipsoid. If the requested number of candidates is not more than the
% dimension + 1, this is done by computing the squared distances of partially
% conditionally rounded float vectors to the float vector in the metric of the
% covariance matrix. Otherwise an approximation is used.
%
% Input arguments
%    L,D   : LtDL-decomposition of the variance-covariance matrix of
%            the float ambiguities (preferably decorrelated)
%    ahat  : float ambiguites (preferably decorrelated)
%    ncands: Requested number of candidates (default = 2)
%    factor: Multiplication factor for the volume of the resulting
%            search ellipsoid (default = 1.5)
%
% Output arguments:
%    Chi2  : Size of the search ellipsoid

% ----------------------------------------------------------------------
% File.....: chistart.m
% Date.....: 19-MAY-1999
% Modified.: 05-MAR-2001, by P. Joosten
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ------------------
% --- Initialize ---
% ------------------

if nargin < 4; ncands = 2  ; end;
if nargin < 5; factor = 1.5; end;

n = max(size(ahat));

% ----------------------------------------------------------------------
% --- Computation depends on the number of candidates to be computed ---
% ----------------------------------------------------------------------

if ncands <= n+1;

  % --------------------------------------------------------
  % --- Computation based on the bootstrapping estimator ---
  % --------------------------------------------------------

  Chi = [];
  iQ  = (L \ diag(1./D)) / L';
  
  for k = n:-1:0;

    afloat = ahat;
    afixed = ahat;
  
    for i = n:-1:1;
        
        dw = 0;
        for j = n:-1:i;
            dw = dw + L(j,i) * (afloat(j) - afixed(j));
        end
        
        afloat(i) = afloat(i) - dw;
        if (i ~= k);
            afixed(i) = round (afloat(i));
        else
            tmp   = round(afloat(i));
            afixed(i) = tmp + sign(afloat(i)-tmp);              
        end
        
    end
    Chi = [Chi (ahat-afixed)' * iQ * (ahat-afixed)];

  end

  % ---------------------------------------------------------------
  % --- Sort the results, and return the appropriate number     ---
  % --- Add an "eps", to make sure there is no boundary problem ---
  % ---------------------------------------------------------------

  Chi  = sort(Chi);
  Chi2 = Chi(ncands) + 1d-6;
   
else

  % -----------------------------------------------------
  % An approximation for the squared norm is computed ---
  % ----------------------------------------------------- 
  
  Vn   = (2/n) * (pi ^ (n/2) / gamma(n/2));
  Chi2 = factor * (ncands / sqrt((prod(D)) * Vn)) ^ (2/n);

end;

% ----------------------------------------------------------------------
% End of routine: chistart
% ----------------------------------------------------------------------
