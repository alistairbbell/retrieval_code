%------------------------------------------------------------------------
% NAME:     FFT_Q.m
%
%           Q definition function for fft_retrieval.m
%
% FORMAT:   Q = FFT_Q( Q )
%
% OUT:      Q    the setting structure for Qpack
% IN:       Q    the setting structure for Qpack
% OPTIONAL: -
%
% TEMPLATE: -
% Q-FIELDS: -
%------------------------------------------------------------------------

% HISTORY: 2001.03.25  Created by Patrick Eriksson
%          2005.07.12  Modified by Soohyun Ka

function Q = FFT_Q( Q )

%=== Report levels
Q.QP_LEVEL   = 2;
Q.ARTS_LEVEL = 0;


%=== Directories

tempdir=pwd;
cd ..
top_dir = pwd;
cd ..
src_dir = pwd;
cd(tempdir);

Q.SRC_DIR        = src_dir;
Q.ARTS           = '/opt/arts/arts1/src/arts';
Q.OUT            = [top_dir,'/Out/out'];
Q.TMP_AREA       = '/tmp';
Q.SPECTRO_DIR    = [top_dir,'/Spectroscopy'];
Q.CALCGRIDS_DIR  = [top_dir,'/CalcGrids'];
Q.SENSOR_DIR     = [top_dir,'/Sensor'];
Q.RETRIEVDEF_DIR = [top_dir,'/RetrievDef'];

Q.APRIORI_VMR = [top_dir,'/RetrievDef/apriori'];
Q.APRIORI_PTZ = [top_dir,'/RetrievDef/apriori.ptz.aa'];

Q.SETUP_VMR = sprintf('"%s/RetrievDef/apriori"',top_dir);
Q.SETUP_PTZ = sprintf('"%s/RetrievDef/apriori.ptz.aa"',top_dir);


%=== Species
Q.RETRIEVAL_TAGS   = '"H2O"';
Q.SPECIES_KGRIDS   = '"retrievalgrid.upperatm.aa"';
Q.SPECIES_COVMATS  =  '"sx.apriori.H2O.v17.aa"';
% Q.SPECIES_COVMATS  =  '"sx.apriori.H2O.Onsala.aa"';

Q.OTHER_TAGS       =  [ ...
      '"H2O-SelfContStandardType",' ...
      '"O2-SelfContStandardType",'...
      '"N2-SelfContStandardType",'...
      '"CO2-SelfContPWR93"'...
      ];

%=== Spectroscopy
Q.LINEFILE         = 'lines_hitran06_jpl01_08.aa';
% Q.LINEFILE         = 'nedoluha1995.aa';

Q.LINESHAPE        = 'Voigt_Drayson';
Q.LINESHAPE_FACTOR = 'quadratic';
Q.LINESHAPE_CUTOFF = -1;

Q.CONTINUA         = 'continua';

%=== Hydrostatic eq.
Q.HSE_IN_ON          = 0;
Q.HSE_RETRIEVAL_ON   = 0;
Q.HSE_PREF           = 95000;
Q.HSE_ZREF           = 0;

%=== RTE
Q.PLATFORM_ALTITUDE  = 16000;
Q.STEPLENGTH_RTE     = 500;

Q.GROUND_ALTITUDE    = 16000;
Q.GROUND_EMISSION    = 0;

Q.REFRACTION_ON      = 0;
Q.EMISSION_ON        = 1;

%=== Calculation grids
Q.P_ABS     = 'p_abs.aa';
Q.F_MONO    = 'opt_fmono_ext.aa';
Q.ZA_PENCIL = 'za_pencil.aa';
Q.F_ORDER   = 1;

%=== Sensor

Q.ANTENNA_ON    = 0;

Q.BACKEND_ON    = 1;
Q.BACKEND_FREQS = 'backend_freqs_100MHz.aa';
Q.BACKEND_FILE  = 'FFT_ACQ1_channelresponse.aa';

Q.DSB_ON        = 1;
Q.DSB_FILE      = 'sbfilter.aa';
Q.DSB_FPRIMARY  = 22.235e9;
Q.DSB_LO        = 20.135e9;

%%% Data reduction

%=== Thermal noise
Q.MEASNOISE_DO     = 2;
Q.MEASNOISE_COVMAT = 'sy.aa';

Q.CALINOISE_DO     = 0;

%=== Retrieval/error quantities beside species
Q.TEMPERATURE_DO       = 2;
Q.TEMPERATURE_KGRID    = 'p_abs.aa';
Q.TEMPERATURE_COVMAT   = 'sx.temperature.aa';

Q.POINTING_DO          = 0;
Q.POINTING_PDF         ='gaussian';
Q.POINTING_STDV        = 0.5;
Q.POINTING_DELTA       = -2;

Q.CONTABS_DO           = 0;

Q.POLYFIT_DO           = 0;
Q.POLYFIT_ORDER        = 0;
Q.POLYFIT_COVMATS      = [];

Q.PPOLYFIT_DO           = 0;

Q.SINEFIT_DO         = 0;
Q.SINEFIT_PERIODS    = [];
Q.SINEFIT_COVMATS    = [];

Q.TB_REFLOADS_DO       = 0;
Q.TB_REFLOADS_NOMINAL  = [30 290];
Q.TB_REFLOADS_STDV     = [5 1];

Q.PROPCAL_DO           = 0;
Q.PROPCAL_TB_REF       = 0;
Q.PROPCAL_STDV         = 0.04;

Q.FREQUENCY_DO		= 0;
Q.FREQUENCY_STDV	= 60e3;
Q.FREQUENCY_DELTA	= 60e3;

Q.PSF_DO     = 0;
Q.PSF_FILE   = 'psf.aa';
Q.PSF_STDV   = 20e3;
Q.PSF_DELTA  = 20e3;

%===Data binning
Q.BINNING_ON           = 0;
Q.BINNING_FILE         =  [top_dir,'/Sensor/binning.aa'];

%=== Retrieval parameters
Q.RETRIEVAL_METHOD     = 'oem';

Q.CLS_SPECIES_POS_ON   = 1;
Q.CLS_NONLIN_ON        = 0;


if Q.CLS_NONLIN_ON
  Q.CLS_GA_START_VALUE      = 1000;
  Q.CLS_GA_FAC_WHEN_OK      = 100;
  Q.CLS_GA_FAC_WHEN_NOT_OK  = 10;
  Q.CLS_STOP                = 0.03;
  Q.CLS_GA_UPPER_LIMIT      = 10000;
  Q.CLS_MAX_ITER            = 10;
end

