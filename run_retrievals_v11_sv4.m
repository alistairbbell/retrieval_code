%%% first setup for retrieval with ARTS/Qpack 2
addpath([pwd '/matlab'])
addpath('/home/miawarac/matlab')
mysql_open;


% Initialize paths for Atmlab:
 run( '/opt/arts/atmlab/atmlab/atmlab_init' );

version=1.1;
sub_version=4;

delta_t=12; %hours
 bw = 80e6;   % Bandwidth in MHz
 deg_poly=2;
    
 %===== save the retrievals in:
 save_dir='/export/data/miawarac/retrieval_files/retrieval_files_v11_sv4';
 unix(['mkdir ' save_dir]);
 
 % start=datenum(2017,11,29);
 start=mysql(sprintf('select max(ROUND(TO_DAYS(time_max)+1+TIME_TO_SEC(time_max)/86400,6)) from MIAWARA_C.level2_header where version between %f and %f and subversion=%i;',version-0.01,version+0.01,sub_version),'mat');
 
 if today == datenum(2019,05,22)
     start = start +10;
 elseif today == datenum(2019,06,13)
     start = start+1;
 end
     
 start=ceil(start*24/delta_t)/(24/delta_t);

 %stop=today;
 stop=mysql('select max(ROUND(TO_DAYS(time)+1+TIME_TO_SEC(time)/86400,6)) from MIAWARA_C.level1_pol1_2504;','mat');
 
    


 
 disp(sprintf('processing time interval %s to %s, ARTS2, %i MHz',datestr(start,31),datestr(stop,31),bw/1e6));
 disp(sprintf('distorted part of spectrum cut away!\n'));
 
 time=start;
 mission=get_mission_from_mysql(time);
 
 for time=start:delta_t/24:stop
  disp(datestr(time))
    %%%%%%%%% Import measurement data %%%%%%%%%
    if mission==11
    [Y, msm] = get_Y_reunion(time,time+delta_t/24,bw);
    else
    [Y, msm] = get_Y(time,time+delta_t/24,bw);
    end


%    [Y, msm] = get_Y_reunion(time,time+delta_t/24,bw); 
    % [Y, msm] = get_Y(time,sigma,bw);    
    if isempty(msm)
    disp('No measurements found')
    %time=time+.1;
    else


    %%%%%%%%% define forward model and retrieval settings %%%%%%%%%
    [Q, O] = setup_Q_O(Y,msm, extra_info);
    
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
    %%% Check that all frequencies are OK
%     if Q.SENSOR_RESPONSE.MIXER_DO
%         Y.F = Q.SENSOR_RESPONSE.F_BACKEND;
%     end
    disp('Checking for consistency...')
    if ~qp2_check_f( Q, Y, 1e3 );
        disp( 'Some mismatch between Q.F_BACKEND and frequencies of spectra.' );
        return
    end

	%%%%%%%%% Perform inversion %%%%%%%%%
	L2 = qpack2(Q,O,Y);

	%%%%%%%%% Add some information and save %%%%%%%%%
    L2.max_time=msm.max_time;
    L2.min_time=msm.min_time;
    L2.tint = msm.tint;
    if isfield(msm,'tau'), 
        L2.tau=nanmean(msm.tau(:,2));
        L2.mirror_elevation=nanmean(msm.mirror_elevation);
        L2.mirror_ele_max=max(msm.mirror_elevation);
        L2.mirror_ele_min=min(msm.mirror_elevation);
    elseif isfield(msm,'tau1')&&~isempty(msm.tau2)&&~isempty(msm.tau1), 
        L2.tau=nanmean([msm.tau1(:,2); msm.tau2(:,2)]);
        L2.mirror_elevation=nanmean([msm.mirror_ele1; msm.mirror_ele2]);
        L2.mirror_ele_max=max([msm.mirror_ele1; msm.mirror_ele2]);
        L2.mirror_ele_min=min([msm.mirror_ele1; msm.mirror_ele2]);
    elseif isfield(msm,'tau1')&&isempty(msm.tau2)&&~isempty(msm.tau1), 
        L2.tau=nanmean(msm.tau1(:,2));
        L2.mirror_elevation=nanmean(msm.mirror_ele1);
        L2.mirror_ele_max=max([msm.mirror_ele1]);
        L2.mirror_ele_min=min([msm.mirror_ele1]);
    elseif isfield(msm,'tau1')&&isempty(msm.tau1)&&~isempty(msm.tau2), 
        L2.tau=nanmean(msm.tau2(:,2));
        L2.mirror_elevation=nanmean(msm.mirror_ele2);
        L2.mirror_ele_max=max([msm.mirror_ele2]);
        L2.mirror_ele_min=min([msm.mirror_ele2]);
    end
  
    L2.Y = Y;
    L2.Q = Q;
    L2.O = O;
    L2.species1_z=interp1(log10(L2.p_grid),L2.z_field,log10(L2.species1_p));
    L2.species1_T=interp1(log10(L2.p_grid),L2.t_field,log10(L2.species1_p));
    if isempty(L2)==0
 
      
       filename=sprintf('%s/retrieval_%i%.02i%.02i%.02i%.02i.mat',save_dir,L2.year,L2.month,L2.day,L2.hour,L2.minute);
       save(filename,'L2');
       
       disp(sprintf('interting version %.1f sv %i',version, sub_version))
       insert_level2_into_mysql(L2,version,sub_version)
       
       %===== the beginning of the next retrieval period is given by the max_time of the current msm
       %time=L2.max_time;
 
    else
 
 
 
    end
    end
 
 end


