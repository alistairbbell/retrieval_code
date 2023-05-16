function save_level2_miawara(obj, filename)

% Saves the matlab object of retrieved fields in a netcdf format

fields = {'year', 'month', 'day', 'hour', 'minute', 'second', 'cost'... 
    'cost_x', 'cost_y', 'species1_p', 'species1_x', 'species1_e'...
    'species1_es', 'species1_mr','species1_z', 'y', 'yf' ...
    'species1_xa', 'converged', 'J', 'f'};

y_fields = {'LONGITUDE', 'LATITUDE' };

array_vals = containers.Map();
array_vals('hour') = [3 2 1];
array_vals('min') =[44 52 66];


disp(fieldnames(obj))
for i = 1:numel(fields)
    fieldname = fields{i};
    temp = double(horzcat(obj(:).(fieldname)));
    disp('size(temp')
    disp(size(temp))
    array_vals(fieldname) = temp;
end

for i = 1:numel(y_fields)
    templist = [];
    fieldname = y_fields{i};
        for j = 1:length(obj)
        templist(end+1) = double(obj(j).Y.(fieldname));
    end
    array_vals(fieldname) = templist;
end

if  exist(filename)
    delete(filename)
    fprintf('Rewriting existing file name: %s', filename)
end

disp('species1_p')
disp(length(obj))
disp(length(obj(1).species1_p))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scientific Dataset 

datetimelist = zeros(length(obj),1);
disp(length(datetimelist))

for i = 1:length(obj)
    dtlist = [obj(i).year obj(i).month obj(i).day obj(i).hour obj(i).minute obj(i).second];
    dt_obs = datetime(dtlist);
    datetimelist(i) = convertTo(dt_obs, 'posixtime');
end

% Coordinates variable (enable 'netcdf4' format)
nccreate(filename,'/time','Dimensions',{'time',length(datetimelist)},'Datatype','double','Format','netcdf4');
nccreate(filename,'/pressure','Dimensions',{'pressure',length(obj(1).species1_p)},'Datatype','double','Format','netcdf4');
nccreate(filename,'/frequency','Dimensions',{'frequency',length(obj(1).f)},'Datatype','double','Format','netcdf4');

%%%%%%%%%%%%%%%%%
% Some variables to help geolocate the file:
% these are scalar variables as they do not vary in time
nccreate(filename,'/lat','Dimensions',{'time',length(datetimelist)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/lon','Dimensions',{'time',length(datetimelist)},'Datatype','single','FillValue',-9999);

% time variables
nccreate(filename,'/z','Dimensions',{'time', length(datetimelist), 'pressure', length(obj(1).species1_p)},'Datatype','single','FillValue',-9999);
%nccreate(filename,'/t_field','Dimensions',{'time', Inf, 'pressure', length(obj.p_grid)},'Datatype','single','FillValue',-9999);

%measurement variables
nccreate(filename,'/y','Dimensions',{'time', length(datetimelist), 'frequency', length(obj(1).f)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/yf','Dimensions',{'time', length(datetimelist), 'frequency', length(obj(1).f)},'Datatype','single','FillValue',-9999);

nccreate(filename,'/q','Dimensions',{'time', length(datetimelist),'pressure', length(obj(1).species1_p)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/q_a','Dimensions',{'time', length(datetimelist),'pressure', length(obj(1).species1_p)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/tau','Dimensions',{'time', length(datetimelist), 'frequency', length(obj(1).f)},'Datatype','single','FillValue',-9999);

%minimisation variables
nccreate(filename,'/cost','Dimensions',{'time', length(datetimelist)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/cost_x','Dimensions',{'time', length(datetimelist)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/cost_y','Dimensions',{'time', length(datetimelist)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/convergence_flag','Dimensions',{'time', length(datetimelist)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/measurement_response','Dimensions',{'time', length(datetimelist),'pressure', length(obj(1).species1_p)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/q_err','Dimensions',{'time', length(datetimelist),'pressure', length(obj(1).species1_p)},'Datatype','single','FillValue',-9999);
nccreate(filename,'/q_err_smooth','Dimensions',{'time', length(datetimelist),'pressure', length(obj(1).species1_p)},'Datatype','single','FillValue',-9999);

%minimisation matrices
nccreate(filename,'/J','Dimensions',{'time', length(datetimelist),...
    'pressure', length(obj(1).species1_p), 'frequency',length(obj(1).f)},'Datatype','single','FillValue',-9999);

%mwr_observation
%nccreate(filename,'/frequency','Dimensions',{'time', Inf, 'frequency', length(obj.f)},'Datatype','single','FillValue',-9999);
%nccreate(filename,'/BT','Dimensions',{'time', Inf, 'frequency', length(obj.f)},'Datatype','single','FillValue',-9999);
%nccreate(filename,'/t_noise','Dimensions', {'time', Inf, 'frequency', length(obj.f)},'Datatype','single','FillValue',-9999);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scientific Dataset (spectrometer1,2,...)
disp('size(datetimelist)')
disp(size(datetimelist))

% Coordinate variables, directly with their attributes
ncwrite(filename,'/time',datetimelist);
ncwriteatt(filename,'/time','units', 's since 1 Jan 1970');
ncwriteatt(filename,'/time','description', ...
    'mean time of beginning of all antenna measurements for this cycle');

ncwrite(filename,'/lat',array_vals('LATITUDE'));
ncwriteatt(filename,'/lat','units', 'degrees North');
ncwriteatt(filename,'/lat', ...
    'description','latitude of the instrument at the time of observation');

ncwrite(filename,'lon',array_vals('LONGITUDE'));
ncwriteatt(filename,'/lon','units', 'degrees North');
ncwriteatt(filename,'/lon','description', ...
    'longitude of the instrument at the time of observation');

%%%%%%%%%%%%%%%%%
% Geolocation variables
%disp('size(species1_p)')
%disp(size(array_vals('species1_p')))

p2d = array_vals('species1_p');
ncwrite(filename,'/pressure', p2d(:,1));
ncwriteatt(filename,'/pressure','units', 'Pa' );
ncwriteatt(filename,'/pressure','description','Pressure grid of retrieved levels');

frequency_vector = array_vals('f');
ncwrite(filename,'/frequency', frequency_vector(:,1));
ncwriteatt(filename,'/frequency','units', 'Hz' );
ncwriteatt(filename,'/frequency','description','Frequency grid of observed spectra');

ncwrite(filename,'/z', array_vals('species1_z')');
ncwriteatt(filename,'/z','units', 'm (asl)' );
ncwriteatt(filename,'/z','description','Height grid of retrieved levels');

ncwrite(filename,'/q', array_vals('species1_x')');
ncwriteatt(filename,'/q','units', 'ppm' );
ncwriteatt(filename,'/q','description','Retrieved water vapour profile');

ncwrite(filename,'/q_err', array_vals('species1_e')');
ncwriteatt(filename,'/q_err','units', 'ppm' );
ncwriteatt(filename,'/q_err','description','Total retrieval (observation + smoothing) error for the water vapour');

ncwrite(filename,'/q_err_smooth', array_vals('species1_es')');
ncwriteatt(filename,'/q_err_smooth','units', 'ppm' );
ncwriteatt(filename,'/q_err_smooth','description','Smoothing error for water vapour profile');

ncwrite(filename,'/measurement_response', array_vals('species1_mr')');
ncwriteatt(filename,'/q','description','Measurement response to change in water vapour content');

ncwrite(filename,'/yf', array_vals('yf')');
ncwriteatt(filename,'/q','units', 'ppm' );
ncwriteatt(filename,'/q','description','Retrieved water vapour profile');

ncwrite(filename,'/y', array_vals('y')');
ncwriteatt(filename,'/q','units', 'ppm' );
ncwriteatt(filename,'/q','description','Retrieved water vapour profile');



disp('size(species1_x)')
disp(size(array_vals('species1_x')))

ncwrite(filename,'/q_a', array_vals('species1_xa')');
ncwriteatt(filename,'/q_a','units', 'ppm' );
ncwriteatt(filename,'/q_a','description','A priori water vapour profile');

ncwrite(filename,'/cost', array_vals('cost'));
ncwriteatt(filename,'/cost','units', '' );
ncwriteatt(filename,'/cost','description','Total cost function value');

ncwrite(filename,'/cost_x', array_vals('cost_x'));
ncwriteatt(filename,'/cost_x','units', '' );
ncwriteatt(filename,'/cost_x','description', ...
    'Cost function value for prior term');

ncwrite(filename,'/cost_y', array_vals('cost_y'));
ncwriteatt(filename,'/cost_y','units', '' );
ncwriteatt(filename,'/cost_y','description', ...
    'Cost function value for measurement term');
disp(['File saved as: ' filename])

ncwrite(filename,'/convergence_flag', array_vals('converged'));
ncwriteatt(filename,'/cost_y','description', ...
    'Convergence of profile in retrieval, -1 = ...');
disp(['File saved as: ' filename])
end
