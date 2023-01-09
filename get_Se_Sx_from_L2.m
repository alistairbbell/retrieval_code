% Derive Se and Sx from L2-structure
% Input: L2       L2-structure
%        full_sx  0: only abs-species Sx, 1: include also polyfit-Sx
function [Se,Sx] = get_Se_Sx_from_L2(L2,full_sx)
if nargin==1
    full_sx = 0;
end

if length( L2.Y.TNOISE )>1
    [i,j,s] = find( L2.Q.TNOISE_C );
    for k = 1 : length(i)
        s(k) = prod( L2.Y.TNOISE([i(k) j(k)]) ) * s(k);
    end
    Se = sparse(i,j,s);    
else
    Se    = (L2.Y.TNOISE^2) .* L2.Q.TNOISE_C;
end
Se = full(Se);

if full_sx&&L2.Q.POLYFIT.RETRIEVE
    n              = length(L2.species1_x);    
    nt             = n+L2.Q.POLYFIT.ORDER+1;
    Sx             = sparse( nt, nt );
    Sx(1:n,1:n)    = L2.Q.ABS_SPECIES(1).SX;
    for i=0:L2.Q.POLYFIT.ORDER
        Sx(n+i+1,n+i+1) = eval(sprintf('L2.Q.POLYFIT.SX%d;',i));
    end

else
    Sx             = L2.Q.ABS_SPECIES(1).SX;
end
Sx = full(Sx);
end