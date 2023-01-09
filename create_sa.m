% function sa = create_sa(corr_fun,cutoff,clength,x0,sx0,xf,x,val)
% Calculate sa-structure used to setup Sx-covariance-matrix in Qpack, 
% which is typically processed in the following way:
% 
%  write_datafile([Q.RETRIEVDEF_DIR,'/sx.apriori.H2O.aa'],sa,'aomatrix');
%
%    where sa is a structure with the following content:
%    sa{1} = [correlation function, correlation cutoff]
%    sa{2} = [altitude-grid, standard deviation, correlation length]
%
%    correlation functions: 0  no correlation, diagonal matrix. The corrlength
%                              and cutoff are of no importance here (see below)
%                           1  linearly decreasing to 0 (tenth function)
%                           2  exponential
%                           3  gaussian   
%
%    correlation cutoff:    correlations below this value are set to 0
%
%    altitude-grid*:        pressure, typically in decades
%                            (log10(pressure)), bottom->top
%    standard deviation     standard deviation of the retrieved quantity at
%                           this altitude
%    correlation length     altitude difference, in which the correlation
%                           goes below 1/e (in the same unit as the alt.-grid!!!))
%
% * If x0/x goes only down to surface, an additional level at 6 decades
%   (p=10^6 Pa) is added at the beginning of the grid to prevent precision
%   problems. Such problems occur because setup_covmatrix (subfunction of
%   sfromfile) checks, if the p-grid in sa covers at least the range of the
%   retrieval grid. If p0(p_grid) is only 1/1000 hPa higher than p0(sa), 
%   this leads to an interrupt of the Sx-calculation.
%
% INPUT:
% corr_fun: correlation function
% cutoff:   correlation cutoff
% clength:  correlation length (either one value, which is then for the 
%           entire grid or a vector in the same size as x0 resp. x)
% x0:       pressure grid(Pa) (= x, if no xf,x defined in input)
% sx0:      std dev.
% xf:       "known" points (fixed points) to fix x_ret   [OPTIONAL]
% x:        final p-grid (retrieval grid)                [OPTIONAL]
% val:      std dev. at fixed points (default: 0.01)
%
%
% OUTPUT: 
% sa        sa-structure as described above
%
% Rene Bleisch 2.2011, rene.bleisch@iap.unibe.ch
%
function sa = create_sa(corr_fun,cutoff,clength,x0,sx0,xf,x,val)

if nargin<8
    val = 0.01;
end
if nargin<7
    x = x0;
end

sa{1}=[corr_fun;cutoff];

sa2 = zeros(length(x),3);

sa2(:,1) = log10(x);

if nargin>5
    sa2(:,2) = interp1(log10(x0),sx0,log10(x));
    
    % add known points
    for i=1:length(xf)
        ind = find(x==xf(i));
        if ~isempty(ind)
            sa2(ind,2) = val;
        else
            ind = find(x>xf(i));
            sa2 = [sa2(ind,:); [log10(xf(i)) val 0]; sa2(ind(end)+1:end,:)];
            x   = [x(ind); xf(i); x(ind(end)+1:end)];
        end
    end
    if length(clength)>1
        sa2(:,3) = interp1(log10(x0),clength,log10(x));
    else
        sa2(:,3) = clength;
    end
else
    sa2(:,2) = sx0;
    sa2(:,3) = clength;
end

% add level below surface to prevent precision problems between 
% p_abs(1) and p(sa)(1) (see *):  
if sa2(1,1)<6
    sa2 = [[6 sa2(1,2:3)]; sa2];
end

sa{2} = sa2;
