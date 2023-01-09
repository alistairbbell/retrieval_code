%%%%%%%%% returns a time interpolated climatology for given time
% at the zimmerwald observatory. The zonal climatology is derived from
% ECMWF global forcasts between the years 2010 and 2014 (inclusive)
%                                                                   %%%%%%%%
function xa = get_apriori_h2o_from_ecmwf_zimmerwald(date, ecmwf_apriori_file)

%disp(date)
MyDayIdx = day(date);
MyMonthIdx = month(date);
MyYear = year(date);
coef_b = 0.6219907;

if exist(ecmwf_apriori_file)
    disp('Found ecmwf file')
    pressure         = ncread(ecmwf_apriori_file,'/P');
    T         = ncread(ecmwf_apriori_file,'/T');
    q         = ncread(ecmwf_apriori_file,'/q');
    Phi_Surf     = ncread(ecmwf_apriori_file,'/Z');
    time = ncread(ecmwf_apriori_file,'/time');
else
    fprintf('ECMWF file %s not found.\n',ecmwf_apriori_file);
end
format long
disp('size(pressure)')

disp(pressure(1,1,1:10))
w = (q./(ones(size(q))-q));
w = w/coef_b;
disp(w(1,1,1:10))

dategrid = ones(3,1);
for i = 1:3
    mytempdate = datenum(datetime(MyYear, MyMonthIdx, 1 ));
    tempdate = datetime(addtodate(mytempdate, i-2, 'month'),'ConvertFrom','datenum');
    dategrid(i) = month(tempdate);
end

DaysInMonth =  days(datetime(year(date), month(date)+1, 1) - datetime(year(date), month(date), 1));
%disp(DaysInMonth)
%disp('size(q)')
%disp(size(q))

month_extract_idx = 1.5+MyDayIdx/DaysInMonth;
time_extract_idx = (hour(date)/24 + minute(date)/(24*60) + second(date)/(24*3600));
%disp('time_extract_idx')
%disp(time_extract_idx)

%extract_array_q = q(:,dategrid, :);
extract_array_q = vertcat(  w(end,dategrid, :), w(:,dategrid, :));
extract_array_q = vertcat(extract_array_q(:,:, :) ,  w(1,dategrid, :));

%extract_array_p = pressure(:,dategrid, :);
extract_array_p = vertcat(pressure(end, dategrid, :), pressure(:,dategrid, :));
extract_array_p = vertcat(extract_array_p(:,:, :) ,  pressure(1,dategrid, :));

disp('shape(extract_array_p)')
disp(size(extract_array_p))

%disp(time)
y_coords = double(time/24)';
y_coords = [y_coords(end) - 1, y_coords];
y_coords = [y_coords, y_coords(2) + 1];

x_coords = double(1:3)';
z_coords = double(1:length(Phi_Surf))';

yq = double(time_extract_idx);
xq = double(month_extract_idx);

%disp('size q')
%disp(extract_array_q)

output_q = interp3(x_coords,y_coords,z_coords, extract_array_q, xq, yq, z_coords);
output_p = interp3(x_coords,y_coords,z_coords, extract_array_p, xq, yq, z_coords);
%disp('x')
%disp(x_coords)
%disp('y')
%disp(y_coords)

xa(:, 1) = output_p;
xa(:, 2) = output_q;
%disp(xa)

end
