% logscale the jacobian and recalculate G and A 
function [Js,Gs,As] = scale_ret(J,x,Se,Sx,dir)
if nargin==4
    dir = 'scale';
elseif nargin<3
    fprintf('Too less arguments!\n');
    fprintf('    [Js,Gs,As] = scale_ret(J,x,Se,Sx[,dir])\n');
    Js   = [];
    Gs   = [];
    As   = [];
    return
end

n        = length(x);

% Calculate Seinv and Sxinv
Seinv    = Se \ speye(size(Se));
Sxinv    = Sx \ speye(size(Sx));

% Scale J
Js = J;

if strcmp(dir,'scale')
    Js(:,1:n) = J(:,1:n)./ repmat( x', size(J,1), 1 );
elseif strcmp(dir,'unscale')
    Js(:,1:n) = J(:,1:n).*repmat( x', size(J,1), 1 );
end

% Recalculate G and A
JtSeinv  = Js' * Seinv;
SJSJ     = Sxinv + JtSeinv * Js;
Gs       = ( SJSJ \ speye(size(Sxinv)) ) * JtSeinv;
As       = Gs*Js;

end