%==========================================================================
% NAME      | get_Y_miawara.m
% TYPE      | Script
% AUTHOR(S) | Alistair Bell
% CREATION  | 01.2023
%           |
% ABSTRACT  | Function to import the l1b measuerement information 
%           | (calibrated and integrated Tb) and assign the Y structure 
%           | through which the retrievals will be passed to the Qpack 
%           | optimal estimation function. 
%           | 
%           | 
% ARGUMENTS | INPUTS: - time
%           |         - extra_info (structure)
%           |
%           | OUTPUTS: Y, msm, extra_info
%           |
% CALLS     | extra_info.get_msm(time_str,extra_info)
%           |
%==========================================================================

function [Y, msm,extra_info ] = get_Y_miawara(time, extra_info)

disp('-> loading measurement')

%definitions for frequency vector
bw_obs = extra_info.bw_obs;
f_cen = extra_info.f_cen;
bw_mn = extra_info.bw_measurement_noise;
fshift = extra_info.use_freq_offset;
twoSigma = extra_info.twoSigma;

%initialize msm and Y
msm.y = [];
Y   = qp2_y;

time_str = sprintf('%s_%s_%s',  num2str(year(time),'%02.f') , num2str(month(time),'%02.f') , num2str(day(time),'%02.f') );

msm = extra_info.get_msm(time_str,extra_info);

%define sigma (measurement noise)
for i=1:length(msm)
    if ~isempty(msm(i).y)==1
      %frequency indices of measurement
      frequency_ind = (find(msm(i).f>=(f_cen+fshift-bw_obs/2) & msm(i).f<=(f_cen+fshift+bw_obs/2)));
      frequency_ind = frequency_ind(~ismember(frequency_ind, extra_info.idx_exclude));
      frequency_ind = sort(frequency_ind);

      frequency_ind_sigma = (find(msm(i).f>=(f_cen+fshift-bw_mn/2) & msm(i).f<=(f_cen+fshift+bw_mn/2)));
      frequency_ind_sigma = frequency_ind_sigma(~ismember(frequency_ind_sigma, extra_info.idx_exclude));
      frequency_ind_sigma = sort(frequency_ind_sigma);

      tb_diff = diff(msm(i).y(frequency_ind_sigma)');
      msm(i).sigma=std(tb_diff);
      disp('TNOISE')
      disp(msm(i).sigma)
    else
       disp('no level1 data found, inversion aborted!')
       msm=[];
       return
    end
end

disp('DONE')
for i = 1:length(msm)
    disp('Start Loop Y init:')
    disp(i)
    %%% set up Y as required by qpack2:
    epochdate  = datetime(1970,1,1);
    Y(i).YEAR       = year(msm(i).time(1)+ datenum(epochdate));
    Y(i).MONTH      = month(msm(i).time(1)+ datenum(epochdate));
    Y(i).DAY        = day(msm(i).time(1)+ datenum(epochdate));
    Y(i).HOUR       = hour(msm(i).time(1)+ datenum(epochdate));
    Y(i).MINUTE     = minute(msm(i).time(1)+ datenum(epochdate));
    Y(i).SECOND     = second(msm(i).time(1)+ datenum(epochdate));
    Y(i).LATITUDE   = msm(i).latitude;
    Y(i).LONGITUDE  = msm(i).longitude;
    
    [time_values, time_idx] = min(abs( msm(i).time - msm(i).meteo_time(:) ) );
    disp(time_idx)
    
    Y(i).HSE_Z      = msm.altitude;
    Y(i).Z_PLATFORM = extra_info.proxy_surf_height; % Platform altitude
    Y(i).ZA  = 0; % Zenith angle of the measurement
    Y(i).F =  msm(i).f; %msm.f-100e3 this field is now obsolete;
    Y(i).F = Y(i).F(frequency_ind);
    Y(i).HSE_P = msm(i).HSE_P(time_idx)*100; %corrected pressure in pa

    %corrected for Tb vector in wrong direction
    Y(i).Y = msm(i).y;
    Y(i).Y = Y(i).Y(frequency_ind);
    
    disp('size(Y(i).Y)')
    disp(size(Y(i).Y))
    
    %%% limit noise to noise >=0.01K!
    %define channel at which noise rises to 2*sigma
    bw_1chan = (max(Y(i).F) - min(Y(i).F)) /(length(Y(i).F));
    %sigma2_chan = int32(extra_info.twoSigma/bw_1chan);
    sigma2_chan = extra_info.twoSigma/bw_1chan;

    if msm(i).sigma > extra_info.sigmaMin
        temp_noise  =  exponential_errors(length(Y(i).Y)/2, sigma2_chan , msm(i).sigma,2)+extra_info.additional_noise;
    else
        %noise is otherwise limited to sigmaMin
        temp_noise  =  exponential_errors(length(Y(i).Y)/2, sigma2_chan ,extra_info.sigmaMin, 2)+extra_info.additional_noise;
    end
    
    %Noise is defined by the function in a symmetric way
    if rem(length(Y(i).Y),2) == 0
      temp_noise2 = [flipud(temp_noise); temp_noise];
    else
      temp_noise2 = [flipud(temp_noise); temp_noise(2:end)];
    end
     
    %Noise temperature assigned to the Y structure based on values obtained
    Y(i).TNOISE = temp_noise2;
    % All measuerments are used - where measuerments are unreliable (index
    % set in extra_info.noise_inf_idx) noise is set to a very large number
    Y(i).TNOISE(extra_info.noise_inf_idx) = 9999999.0;
    disp('length(Y.TNOISE)')
    disp(size(Y(i).TNOISE))
end

function output = exponential_errors(N_fchans, chn_double, err,base)
    %Nchans - number of channels errors output for
    %N_double - channel number at which error is doubled
    %sigma - error at first channel
    %returns an array of exponentially increasing values
    output = zeros(N_fchans,1);
    doub_fac = log(2)/log(base);
    for j =1:N_fchans
        output(j) = base^(j*doub_fac/chn_double);
    end
    output = output*err;
end

end
