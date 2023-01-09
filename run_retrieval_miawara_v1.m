%%% first setup for retrieval with ARTS/Qpack 2.4
%addpath(['/home/alistair/MIAWARA_ret/old_miawara_calibration/matlab'])
%addpath('//home/alistair/')
%addpath('/home/alistair/ARTS/testing/atmlab/atmlab-2.2.22/retrieval/')
% Initialize paths for Atmlab:
addpath('/home/alistair/ARTS/testing/atmlab-2.4.2/')
run( '/home/alistair/ARTS/testing/atmlab-2.4.2/atmlab/atmlab_init.m' )

version=1;
sub_version=4;

delta_t=24;  %hours
bw = 80e6;   % Bandwidth in MHz
deg_poly=2;

%===== information needed for retrievals=====================:
extra_info.proxy_surf_height = 1000; %minimum retrieval height in m
extra_info.l1_savename = 'MIAWARA_level1a_AC240_' ;
extra_info.l1_file_ext = '.nc' ;
extra_info.save_dir='//home/alistair/export/data/miawara/Retrievals_ppmv_conv' ;
extra_info.l1_dir = '/home/alistair/export/data/miawara/MIA_calibration_GROSOM-harmo/' ;
extra_info.top_dir = '/home/alistair/MIAWARA_ret/old_miawara_calibration/retrieval_data' ;
extra_info.ecmwf_dir = '/storage/tub/atmosphere/ecmwf/locations/Bern/' ;
extra_info.ecmwf_filename = 'ecmwf_oper_v2_BERN_' ;
extra_info.ecmwf_coeff = '/home/alistair/MIAWARA_ret/extra_files/table_ECMWF.csv';
extra_info.error_list = strings(0);

%monthly a priori defined here
extra_info.a_priori_filename = '/home/alistair/MIAWARA_ret/a_priori/ecmwf_2010_2015_3_9_15_21h.nc' ;

unix(['mkdir ' extra_info.save_dir]);
start=datenum(2022,10, 20 );
stop=datenum(2022,10, 20 );

start = ceil(start*24/delta_t)/(24/delta_t);

disp('start')
disp(start)
stop = ceil(stop*24/delta_t)/(24/delta_t);
disp('stop')
disp(stop)

disp(sprintf('processing time interval %s to %s, ARTS2, %i MHz',datestr(start,31),datestr(stop,31),bw/1e6));
disp(sprintf('distorted part of spectrum cut away!\n'));

time=start;
 for time=start:delta_t/24:stop
  %disp(datestr(time))
    %%%%%%%%% Import measurement data %%%%%%%%
    try
        [Y, msm] = get_Y_miawara(time,bw, extra_info);
        disp('msm fieldnames')
        disp(fieldnames(msm))
        disp(msm(1).time)
        if isempty(msm)
        disp('No measurements found')
        %time=time+.1;
        else
    
    
        for i = 1:length(msm)
            %%%%%%%%% define forward model and retrieval settings %%%%%%%%%
            disp('Profile number:')
            disp(i)
            [Q, O] = setup_Q_O(Y(i),msm(i), extra_info);
            
            if isempty(Q)
                %time=time+.1;
                continue
            end
            Q.POLYFIT.RETRIEVE        = true;
            Q.POLYFIT.ORDER           = deg_poly;
            Q.POLYFIT.L2              = true;
            Q.POLYFIT.SX0             = 1^2; 
            Q.POLYFIT.SX1             = 0.5^2; 
            Q.POLYFIT.SX2             = 0.2^2;
            Q.POLYFIT.SX3             = [];%0.1^2;
            Q.POLYFIT.SX4             = [];%0.1^2; 
            Q.POLYFIT.SX5             = [];%0.1^2; 
            Q.POLYFIT.SX6             = [];%0.1^2;
            Q.MBLOCK_DLOS_GRID    = 0;
            
            %%% Check that all frequencies are OK
        
            Q.SENSOR_RESPONSE.MIXER_DO = cellfun(@(x) ~isempty(x) && x~=0,Q.SENSOR_RESPONSE.MIXER_DO);
            disp('Q.Z.ATMDATA')
            disp(Q.Z.ATMDATA)
            disp(Q.SENSOR_RESPONSE.MIXER_DO)
            disp('Q.SENSOR_RESPONSE.MIXER_DO')
        
            %if Q.SENSOR_RESPONSE.MIXER_DO
            Y(i).F = Q.SENSOR_RESPONSE.F_BACKEND;
            %end
            disp('Checking for consistency...')
        
            disp( 'Q.SENSOR_RESPONSE.F_BACKEND')
            disp( size(Q.SENSOR_RESPONSE.F_BACKEND))
            
            if ~qp2_check_f( Q, Y(i), 1e3 );
                disp( 'Some mismatch between Q.F_BACKEND and frequencies of spectra.' );
                return
        end
        disp('====START INVERSION======')
        try
	        %%%%%%%%% Perform inversion %%%%%%%%%
	        L2 = qpack2(Q,O,Y(i));
            %L2 = qpack2(Q,O,Y(i));

	        %%%%%%%%% Add some information and save %%%%%%%%%
            L2.max_time=msm(i).max_time;
            L2.min_time=msm(i).min_time;
            L2.tint = msm(i).tint;
            if isfield(msm,'tau') 
               % L2.tau=mean(msm.tau(:,2),'omitnan');
                L2.tau=mean(msm(i).tau(:),'omitnan');
                L2.mirror_elevation= mean(msm(i).mirror_elevation,'omitnan');
                L2.mirror_ele_max=max(msm(i).mirror_elevation);
                L2.mirror_ele_min=min(msm(i).mirror_elevation);
            elseif isfield(msm(i),'tau1')&&~isempty(msm(i).tau2)&&~isempty(msm(i).tau1), 
                L2.tau=mean([msm(i).tau1(:,2); msm(i).tau2(:,2)],'omitnan' );
                L2.mirror_elevation=nanmean([msm(i).mirror_ele1; msm(i).mirror_ele2]);
                L2.mirror_ele_max=max([msm(i).mirror_ele1; msm(i).mirror_ele2]);
                L2.mirror_ele_min=min([msm(i).mirror_ele1; msm(i).mirror_ele2]);
            elseif isfield(msm(i),'tau1')&&isempty(msm(i).tau2)&&~isempty(msm(i).tau1), 
                L2.tau=mean(msm(i).tau1(:,2), 'omitnan');
                L2.mirror_elevation=mean(msm(i).mirror_ele1, 'omitnan');
                L2.mirror_ele_max=max([msm(i).mirror_ele1]);
                L2.mirror_ele_min=min([msm(i).mirror_ele1]);
            elseif isfield(msm(i),'tau1')&&isempty(msm(i).tau1)&&~isempty(msm(i).tau2), 
                L2.tau=mean(msm(i).tau2(:,2), 'omitnan');
                L2.mirror_elevation=mean(msm(i).mirror_ele2, 'omitnan');
                L2.mirror_ele_max=max([msm(i).mirror_ele2]);
                L2.mirror_ele_min=min([msm(i).mirror_ele2]);
            end
            
            L2.Y = Y(i);
            L2.Q = Q;
            L2.O = O;
            L2.species1_z=interp1(log10(L2.p_grid),L2.z_field,log10(L2.species1_p));
            L2.species1_T=interp1(log10(L2.p_grid),L2.t_field,log10(L2.species1_p));
            if isempty(L2)==0
               
               filename=sprintf('%s/retrieval_%i%.02i%.02i%.02i%.02i.mat',extra_info.save_dir,L2.year,L2.month,L2.day,L2.hour,L2.minute);
               disp(filename)
               save(filename,'L2');
               
               
               %===== the beginning of the next retrieval period is given by the max_time of the current msm
               %time=L2.max_time;
         
            else
                disp('L2 is empty')
         
            end
    catch ME
        extra_info.error_list(1) = strcat(ME.message,string(datetime((msm(i).time + datenum(1970,1,1)),'ConvertFrom','datenum')))  ;
    end
    disp(extra_info.error_list)
      end
        end
    end
end


