function mission=get_mission_from_mysql(time)
mysql_open;
mysql('use MIAWARA_C;');
cmd=sprintf('select mission_ID from mission where start_date<"%s" and end_date>"%s"',datestr(time+1,31),datestr(time-1,31));
mission=mysql(cmd,'mat');
