% function msm = loadmsm(M)
%
% This function loads the measurments of MIAWARA (level1 data)
% from MySQL, according to the information in the M structure.
% The M-structure contains information such as start time, stop time,
% spectrometer, type of integration, etc...
%
% Example for an M-Structure:
%
%   M.spectrometer = 'FFT';               % Which spectrometer to use ('FFT', 'CTS' or 'BEAM')
%   M.bandwidth = 80e6;                   % Bandwidth in [Hz]
%   M.int_type = 'fixed_sigma_low_noise'; % Integration type:
%                                           - 'fixed_sigma'
%                                           - 'fixed_sigma_low_noise'
%                                           - 'fixed_period'
%                                           - 'fixex_period_low_noise'
%                                           - 'from_file'
%   M.starttime = '2011-01-02 00:30:00';  % Start time of the measurements
%                                           ('yyyy-mm-dd HH:MM:SS')
%   M.stoptime = '';                      % Stop time of the measurements,
%                                           optional depending on integration type
%                                           ('yyyy-mm-dd HH:MM:SS')
%   M.sigma0 = 0.01;                      % Target noise level [K],
%                                           optional depending on integration type
%   M.hour_interval = [0 3];              % Hour interval, only measurements within
%                                           the given hours are integrated, optional
%   M.filename = 'miawara_msm.mat'        % File name of .mat-File to load if integration
%                                           type if 'from_file'.
%
% DS, 2012-01-09
%
function msm = loadmsm(M)

if isfield(M,'hour_interval') && ~isempty(M.hour_interval) && ~strcmp(M.int_type,'fixed_period')
    error('If hour_interval is set, integration type must be fixed_period!')
elseif ~isfield(M,'hour_interval')
    M.hour_interval = [];
end

switch M.int_type
    case 'fixed_sigma'
        msm = integrate_miawara_level1_fixed_sigma(M.starttime,M.sigma0,M.spectrometer);
    case 'fixed_period'
        msm = integrate_miawara_level1_fixed_period(M.starttime,M.stoptime,M.spectrometer,M.ref_angle,M.hour_interval);
    case 'fixed_sigma_low_noise'
        msm = integrate_miawara_level1_fixed_sigma_lownoise(M.starttime,M.sigma0,M.spectrometer);
    case 'fixed_period_low_noise'
        msm = integrate_miawara_level1_fixed_period_lownoise(M.starttime,M.stoptime,M.spectrometer);
    case 'fixed_period_tropo'
        msm = get_miawara_level1_fft_total_power(M.starttime,M.stoptime,1,M.time);
    case 'from_file'
        if isfield(M,'filename')
            load(M.filename)
        else
            error('No file name given!');
        end
    otherwise
        error('Type of integration not recognised!')
end

if ~isempty(msm) && ~isempty(msm.y)
    if strcmp(M.spectrometer,'CTS')
        % Remove the bad channels of the CTS spectrometer
        msm.y = msm.y(cts_good_channels);
        msm.f = msm.f(cts_good_channels);
    end
    if strcmp(M.spectrometer,'BEAM')
        % Remove the center channel of the BEAM spectrometer
        msm.y(1025) = [];
        msm.f(1025) = [];
    end
    
    f_ind = flipud(find(msm.f>=22.23508e9-M.bandwidth/2 & msm.f<=22.23508e9+M.bandwidth/2));
    msm.y = msm.y(f_ind);
    msm.f = msm.f(f_ind);
else
    msm = [];
end
