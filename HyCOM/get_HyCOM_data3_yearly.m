% Download HyCOM reanlysis data using OPeNDAP
% I am downloading data from the ftp site (see /projects/rsmas/kirtman/rxb826/DATA/HYCOM_Reanlysis/sfccur/get_data.ksh) however it's painfully slow.
% Test using OPenDAP to get the variables and levels I want without having to download the entire file.
% Using the data in /projects/rsmas/kirtman/rxb826/DATA/HYCOM_Reanlysis/sfccur/ as a NetCDF template for which to write the data to
%
% Exp 19.0 has data from 19940101 - 19950731
% OPeNDAP website is http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0/3hrly
% The first time point (0) is 19921002_t000 and the next time point (1) is
% 19921002_t003
% For 19940101 time starts at ...
%
% Exp 19.1 has data from 19950801 - 20121231
%
% Don't forget 'module load matlab'
%
% If freezes for more than 30 minutes kill and start again
%
% Jan 10th 2017
clear variables; close all
%profile on
cd '/Volumes/SAMSUNG/WORK/POSTDOC_RSMAS_2016/MATLAB/get_HyCOM_data'

opendapyear=2008;
% Change these
%startdate=19950419;
startdate=str2double([num2str(opendapyear),'1012']);
%enddate=19950801;
enddate=str2double([num2str(opendapyear),'1231']);

%exp='19.0';expnum='190';
exp='19.1';expnum='191'; % For data after 19950801
hour='15'; % 00, 03, 06, 09. 12, 15, 18 or 21
hournum=str2double(hour)/24;

tsdver=0; % 0, 1, 2 or 3
if tsdver==0
    data_url=['http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_',exp,'/',num2str(opendapyear),'/3hrly'];
end
% if tsdver==1
%     data_url=['http://tds1.hycom.org/thredds/dodsC/GLBu0.08/expt_',exp,'/3hrly'];
% end
% if tsdver==2
%     data_url=['http://tds2.hycom.org/thredds/dodsC/GLBu0.08/expt_',exp,'/3hrly'];
% end
% if tsdver==3
%     data_url=['http://tds3.hycom.org/thredds/dodsC/GLBu0.08/expt_',exp,'/3hrly'];
% end

% For when macfusion is not working don't check the file and put it onto
% the harddrive to scp it across later
macfusion=0;
if macfusion==0
    odir='/Volumes/Seagate Backup Plus Drive/Seagate Dashboard/WORK/POSTDOC_RSMAS_2016/DATA/HYCOM_Reanlysis/sfccur/';
else
    odir='/Volumes/Pegasus/HYCOM_Reanlysis/sfccur/';
end

% This is the time the netcdf files have
time_origin='01-Jan-2000 00:00:00';
time_origin_num=datenum(time_origin);

% Belive startdatestr is the same for each year but may need to check
if strcmp(exp,'19.0')
    datastartdatestr='02-Oct-1992 00:00:00';
else
    datastartdatestr=['01-Jan-',num2str(opendapyear),' 00:00:00'];
end
datastartdatenum=datenum(datastartdatestr);

% Field settings
depth1='0';
depth2='4';
depth_sel=['[',depth1,':1:',depth2,']'];
lat1='0';
lat2='2000';
lat_sel=['[',lat1,':1:',lat2,']'];
lon1='0';
lon2='4499';
lon_sel=['[',lon1,':1:',lon2,']'];

% Time points may vary year to year
if strcmp(num2str(opendapyear),'2008')
    ntimepoints='2923';
end    
url_file = [data_url,'?depth[0],lat[1000],lon[2250],time[0:1:',ntimepoints,'],water_u[0:1:',ntimepoints,'][0][1000][2250]'];    

all_time = ncread(url_file,'time');  % 'hours since 2000-01-01 00:00:00'
% Convert all_time to date
all_time_days = all_time / 24;
all_time_datenum = time_origin_num + all_time_days;
all_time_datestr = datestr(all_time_datenum);

% Reference file for which to save netcdf file like
reffile='hycom_uv_GLBu0.08_190_1994010100_t006_sfc10m.nc';

reffiletimeorgin='01-Jan-1850 00:00:00';
reffiletimeorgin_num=datenum(reffiletimeorgin);
% Netcdf to save into reference
finfo=ncinfo([odir,reffile]);
lonfile = ncread([odir,reffile],'lon');
latfile = ncread([odir,reffile],'lat');
Ufile = ncread([odir,reffile],'water_u');

% Loop from startdate to enddate

% Convert these to date strings in format DD-MMM-YYYY (i.d. 1
% https://www.mathworks.com/help/matlab/ref/datestr.html)
tmp=num2str(startdate);
if strcmp(tmp(5:6),'01')
    monstr='Jan';
elseif strcmp(tmp(5:6),'02')
    monstr='Feb';
elseif strcmp(tmp(5:6),'03')
    monstr='Mar';
elseif strcmp(tmp(5:6),'04')
    monstr='Apr';
elseif strcmp(tmp(5:6),'05')
    monstr='May';
elseif strcmp(tmp(5:6),'06')
    monstr='Jun';
elseif strcmp(tmp(5:6),'07')
    monstr='Jul';
elseif strcmp(tmp(5:6),'08')
    monstr='Aug';
elseif strcmp(tmp(5:6),'09')
    monstr='Sep';
elseif strcmp(tmp(5:6),'10')
    monstr='Oct';
elseif strcmp(tmp(5:6),'11')
    monstr='Nov';
elseif strcmp(tmp(5:6),'12')
    monstr='Dec';
end
startdatestr=[tmp(7:8),'-',monstr,'-',tmp(1:4)];
tmp=num2str(enddate);
if strcmp(tmp(5:6),'01')
    monstr='Jan';
elseif strcmp(tmp(5:6),'02')
    monstr='Feb';
elseif strcmp(tmp(5:6),'03')
    monstr='Mar';
elseif strcmp(tmp(5:6),'04')
    monstr='Apr';
elseif strcmp(tmp(5:6),'05')
    monstr='May';
elseif strcmp(tmp(5:6),'06')
    monstr='Jun';
elseif strcmp(tmp(5:6),'07')
    monstr='Jul';
elseif strcmp(tmp(5:6),'08')
    monstr='Aug';
elseif strcmp(tmp(5:6),'09')
    monstr='Sep';
elseif strcmp(tmp(5:6),'10')
    monstr='Oct';
elseif strcmp(tmp(5:6),'11')
    monstr='Nov';
elseif strcmp(tmp(5:6),'12')
    monstr='Dec';
end
enddatestr=[tmp(7:8),'-',monstr,'-',tmp(1:4)];
NumDays=daysact(startdatestr,enddatestr);

filedatenum=datenum(startdatestr)+hournum;

missingdatacounter = 1;
% Creat empty string array for missing data
[row,col] = size(all_time_datestr);
missingdata=char.empty(0,col);

for i = 1:NumDays
%for i = 1:2
    
    % Find filedatanum in all_time_datenum
    timeindex = find(all_time_datenum==filedatenum);
    
    % Check if the netcdf file exits in odir
    % Convert filedatanum into YYYYMMDD00_t0HH
    tmp=datestr(filedatenum,'yyyy-mm-dd HH:MM:SS');
    ncfilename = ['hycom_uv_GLBu0.08_',expnum,'_',tmp(1:4),tmp(6:7),tmp(9:10),'00_t0',tmp(12:13),'_sfc10m.nc'];
    disp(ncfilename)
    if exist([odir,ncfilename], 'file') == 2
        filedatenum = filedatenum + 1; % Increase by one day
    else
        % Check if value is empty (i.e. data is missing)
        if isempty(timeindex) == 1
            % Missing data
            disp(['missing data for ',datestr(filedatenum)])
            missingdata(missingdatacounter,:) = datestr(filedatenum);
            missingdatacounter = missingdatacounter + 1;
            filedatenum = filedatenum + 1; % Increase by one day
        else
            clear all_time url_file U V
            url_file = [data_url,'?depth',depth_sel,',lat',lat_sel,',lon',lon_sel,',','time[',num2str(timeindex),'],water_u[',num2str(timeindex),']',depth_sel,lat_sel,lon_sel,',water_v[',num2str(timeindex),']',depth_sel,lat_sel,lon_sel];
            % ncdisp(url_file);
            U = ncread(url_file,'water_u');
            % Average all the depths of U
            Uavg = mean(U,3);
            
            % Check data looks ok
            
            V = ncread(url_file,'water_v');
            % Average all the depths of U
            Vavg = mean(V,3);
            
            % Work out time as days since 1850-01-01 00:00:00
            filetime = filedatenum - reffiletimeorgin_num;
            
            % Write netcdf file
            fout=ncfilename;
            delete([odir,fout])
            ncwriteschema([odir,fout],finfo);
            ncwrite([odir,fout],'lon',lonfile);
            ncwrite([odir,fout],'lat',latfile);
            ncwrite([odir,fout],'time',filetime);
            ncwrite([odir,fout],'water_u',Uavg);
            ncwrite([odir,fout],'water_v',Vavg);            

            filedatenum = filedatenum + 1; % Increase by one day
        end
    end
    %profile viewer
end

% Save missingdata as text file
% Check if these missing data are on the ftp website
savestr = ['missing_files_',num2str(startdate),'-',num2str(enddate),'_t0',hour,'.txt'];
fileID = fopen(savestr,'w');
formatSpec='%s\n';
[nrows,ncols] = size(missingdata);
for row = 1:nrows
    fprintf(fileID,formatSpec,missingdata(row,:));
end
