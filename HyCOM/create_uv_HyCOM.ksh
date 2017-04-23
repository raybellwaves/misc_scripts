#!/bin/ksh
#
# Interpolate HyCOM data on 0.5x0.5o grid using NCL
# module load netcdf/4.2.1.1 ncl
# bsub < create_uv_HyCOM_submit.sh or ./create_uv_HyCOM.ksh
# to run in command line ncl batch_files/file_0.ncl
# Note can't run other scripts in this directory at this time
#
# 2/22/17

if [ -d 'batch_files' ];then
   :
else
   mkdir batch_files
fi

ddir=/projects/rsmas/kirtman/rxb826/DATA/HYCOM_Reanlysis/sfccur/
odir=/projects/rsmas/kirtman/rxb826/DATA/HYCOM_Reanlysis/sfccur_0.5x0.5/

# counter
fileref=0

for year in 1994 ; do
#for year in 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012; do
   #for month in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
   for month in 01 ; do
      #for day in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 ; do
      for day in 01 ; do
         #for hour in 00 03 06 09 12 15 18 21 ; do
         for hour in 00 ; do

            # Check if in file exists
            fileexists=0
            if [[ -f ${ddir}hycom_uv_GLBu0.08_190_${year}${month}${day}00_t0${hour}_sfc10m.nc ]]; then
               expnum=190
               inputfile=${ddir}hycom_uv_GLBu0.08_${expnum}_${year}${month}${day}00_t0${hour}_sfc10m.nc
               fileexists=1
            fi
            if [[ -f ${ddir}hycom_uv_GLBu0.08_191_${year}${month}${day}00_t0${hour}_sfc10m.nc ]]; then
               expnum=191
               inputfile=${ddir}hycom_uv_GLBu0.08_${expnum}_${year}${month}${day}00_t0${hour}_sfc10m.nc
               fileexists=1
            fi
            if [[ -f ${ddir}hycom_uv_GLBu0.08_19m_${year}${month}${day}00_t0${hour}_sfc10m.nc ]]; then
               expnum=19m
               inputfile=${ddir}hycom_uv_GLBu0.08_${expnum}_${year}${month}${day}00_t0${hour}_sfc10m.nc
               fileexists=1
            fi
            if [ $fileexists -eq 0 ];then
               echo "file ${ddir}hycom_uv_GLBu0.08_19*_${year}${month}${day}00_t0${hour}_sfc10m.nc is missing"
               #exit 0
            else
               rm -rf batch_files/file_${fileref}.ncl
               echo 'load "gsn_code.ncl"' > batch_files/file_${fileref}.ncl
               chmod u+x batch_files/file_${fileref}.ncl
               # Load functions
               echo 'load "gsn_csm.ncl"' >> batch_files/file_${fileref}.ncl
               echo 'load "contributed.ncl"' >> batch_files/file_${fileref}.ncl
               # Setup directories
               echo 'indir = "'${ddir}'"' >> batch_files/file_${fileref}.ncl
               echo 'outdir = "'${odir}'"' >> batch_files/file_${fileref}.ncl  
               echo 'xx1 = new((/1,2001,4500/),float)' >> batch_files/file_${fileref}.ncl
               echo 'xxx1 = new((/1,323,720/),float)' >> batch_files/file_${fileref}.ncl
               echo 'year = "'${year}'"' >> batch_files/file_${fileref}.ncl
               echo 'month = "'${month}'"' >> batch_files/file_${fileref}.ncl 
               echo 'day = "'${day}'"' >> batch_files/file_${fileref}.ncl
               echo 'hour = "'${hour}'"' >> batch_files/file_${fileref}.ncl
               echo 'expnum = "'${expnum}'"' >> batch_files/file_${fileref}.ncl
               # Read in fileG
               echo 'filin = "'${inputfile}'"' >> batch_files/file_${fileref}.ncl
               echo 'filein = addfile(filin,"''r")' >> batch_files/file_${fileref}.ncl
               # U
               echo 'filo = "hycom_u_GLBu_"+expnum+"_"+year+month+day+"00_t0"+hour+"_sfc10m.nc"' >> batch_files/file_${fileref}.ncl
               # Remove the file if exists
               echo 'system ("/bin/rm -f "+outdir+filo)' >> batch_files/file_${fileref}.ncl
               echo 'fout = addfile (outdir+filo , "c")' >> batch_files/file_${fileref}.ncl
               echo 'delete(xx1)' >> batch_files/file_${fileref}.ncl
               # Original resolution (time, lat, lon)
               echo 'xx1 = new((/1,2001,4500/),float)' >> batch_files/file_${fileref}.ncl
               echo 'x1 = filein->water_u' >> batch_files/file_${fileref}.ncl
               # Extract a missing value from this
               echo 'hycom_nanflag = x1(0,0,0)' >> batch_files/file_${fileref}.ncl
               # Put this data into the xx variable
               echo 'xx1 = x1' >> batch_files/file_${fileref}.ncl
               # Add lat
               echo 'xx1!1 = "lat"' >> batch_files/file_${fileref}.ncl
               echo 'xx1&lat = x1&lat' >> batch_files/file_${fileref}.ncl
               # Add lon
               echo 'xx1!2 = "lon"' >> batch_files/file_${fileref}.ncl
               echo 'xx1&lon = x1&lon' >> batch_files/file_${fileref}.ncl
               # Add name of variable
               echo 'xx1@standard_name = "eastward_sea_water_velocity"' >> batch_files/file_${fileref}.ncl               
               echo 'xx1@long_name = "Eastward Water Velocity"' >> batch_files/file_${fileref}.ncl              
               echo 'xx1@units = "m/s"' >> batch_files/file_${fileref}.ncl
               ## Print dimensions of xx1
               ##echo 'print(dimsizes(xx1))' >> batch_files/file_${fileref}.ncl
               # 0.5x0.5 re-gridding
               echo 'delete(xxx1)' >> batch_files/file_${fileref}.ncl 
               echo 'xxx1 = new((/1,323,720/),float)' >> batch_files/file_${fileref}.ncl
               echo 'xi = xx1&lon' >> batch_files/file_${fileref}.ncl
               echo 'yi = xx1&lat' >> batch_files/file_${fileref}.ncl
               echo 'xo = fspan(-180., 179.5, 720)' >> batch_files/file_${fileref}.ncl
               echo 'yo = fspan(-80., 80.,323)' >> batch_files/file_${fileref}.ncl
               echo 'xxx1 = linint2_Wrap (xi, yi, xx1, False, xo, yo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 80oS
               echo 'xxx1(:,0,:) = linint1_Wrap (xi, xx1(:,0,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 79.5oS
               echo 'xxx1(:,1,:) = linint1_Wrap (xi, xx1(:,0,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 79oN
               echo 'xxx1(:,321,:) = linint1_Wrap (xi, xx1(:,2000,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 90oN
               echo 'xxx1(:,322,:) = linint1_Wrap (xi, xx1(:,2000,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 359oE
               echo 'xxx1(:,:,719) = linint1_Wrap (yi, xx1(:,:,4449), False, yo, 0)' >> batch_files/file_${fileref}.ncl
               # Re-do 80oS
               echo 'xxx1(:,0,:) = linint1_Wrap (xi, xx1(:,0,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Re-do 80oN
               echo 'xxx1(:,322,:) = linint1_Wrap (xi, xx1(:,2000,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1 = new((/1,361,720/),float)' >> batch_files/file_${fileref}.ncl
               echo 'yoo = fspan(-90., 90.,361)' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1(:,19:341,:) = xxx1' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1(:,0:18,:) = hycom_nanflag' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1(:,341:360,:) = hycom_nanflag' >> batch_files/file_${fileref}.ncl
              # Add time
               echo 'time = x1&time' >> batch_files/file_${fileref}.ncl
               # Write output file
               echo 'xxxx1!0 = "time"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time = time' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@standard_name = "time"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@long_name_name = "time"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@units = "days since 1850-01-01 00:00:00"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@calendar = "standard"'  >> batch_files/file_${fileref}.ncl
               echo 'xxxx1!1 = "lat"' >> batch_files/file_${fileref}.ncl
               # Add latitude
               echo 'xxxx1&lat = yoo' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@data_type = "float"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@units = "degrees_north"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@axis = "Y"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@standard_name = "latitude"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@valid_min = -90.0' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@valid_max =  90.0' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1!2 = "lon"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@data_type = "float"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@units = "degrees"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@axis = "X"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@standard_name = "longitude"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@valid_min = -180.' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@valid_max = 179.5' >> batch_files/file_${fileref}.ncl
               echo 'fout->water_u = xxxx1' >> batch_files/file_${fileref}.ncl
               # V
               echo 'filo = "hycom_v_GLBu_"+expnum+"_"+year+month+day+"00_t0"+hour+"_sfc10m.nc"' >> batch_files/file_${fileref}.ncl
               # Remove the file if exists
               echo 'system ("/bin/rm -f "+outdir+filo)' >> batch_files/file_${fileref}.ncl
               echo 'fout = addfile (outdir+filo , "c")' >> batch_files/file_${fileref}.ncl
               echo 'delete(xx1)' >> batch_files/file_${fileref}.ncl
               # Original resolution (time, lat, lon)
               echo 'xx1 = new((/1,2001,4500/),float)' >> batch_files/file_${fileref}.ncl
               echo 'x1 = filein->water_v' >> batch_files/file_${fileref}.ncl
               # Extract a missing value from this
               echo 'hycom_nanflag = x1(0,0,0)' >> batch_files/file_${fileref}.ncl
               # Put this data into the xx variable
               echo 'xx1 = x1' >> batch_files/file_${fileref}.ncl
               # Add lat
               echo 'xx1!1 = "lat"' >> batch_files/file_${fileref}.ncl
               echo 'xx1&lat = x1&lat' >> batch_files/file_${fileref}.ncl
               # Add lon
               echo 'xx1!2 = "lon"' >> batch_files/file_${fileref}.ncl
               echo 'xx1&lon = x1&lon' >> batch_files/file_${fileref}.ncl
               # Add name of variable
               echo 'xx1@standard_name = "northward_sea_water_velocity"' >> batch_files/file_${fileref}.ncl               
               echo 'xx1@long_name = "Northward Water Velocity"' >> batch_files/file_${fileref}.ncl              
               echo 'xx1@units = "m/s"' >> batch_files/file_${fileref}.ncl
               ## Print dimensions of xx1
               ##echo 'print(dimsizes(xx1))' >> batch_files/file_${fileref}.ncl
               # 0.5x0.5 re-gridding
               echo 'delete(xxx1)' >> batch_files/file_${fileref}.ncl 
               echo 'xxx1 = new((/1,323,720/),float)' >> batch_files/file_${fileref}.ncl
               echo 'xi = xx1&lon' >> batch_files/file_${fileref}.ncl
               echo 'yi = xx1&lat' >> batch_files/file_${fileref}.ncl
               echo 'xo = fspan(-180., 179.5, 720)' >> batch_files/file_${fileref}.ncl
               echo 'yo = fspan(-80., 80.,323)' >> batch_files/file_${fileref}.ncl
               echo 'xxx1 = linint2_Wrap (xi, yi, xx1, False, xo, yo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 80oS
               echo 'xxx1(:,0,:) = linint1_Wrap (xi, xx1(:,0,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 79.5oS
               echo 'xxx1(:,1,:) = linint1_Wrap (xi, xx1(:,0,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 79oN
               echo 'xxx1(:,321,:) = linint1_Wrap (xi, xx1(:,2000,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 90oN
               echo 'xxx1(:,322,:) = linint1_Wrap (xi, xx1(:,2000,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Replace missing values at 359oE
               echo 'xxx1(:,:,719) = linint1_Wrap (yi, xx1(:,:,4449), False, yo, 0)' >> batch_files/file_${fileref}.ncl
               # Re-do 80oS
               echo 'xxx1(:,0,:) = linint1_Wrap (xi, xx1(:,0,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               # Re-do 80oN
               echo 'xxx1(:,322,:) = linint1_Wrap (xi, xx1(:,2000,:), False, xo, 0)' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1 = new((/1,361,720/),float)' >> batch_files/file_${fileref}.ncl
               echo 'yoo = fspan(-90., 90.,361)' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1(:,19:341,:) = xxx1' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1(:,0:18,:) = hycom_nanflag' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1(:,341:360,:) = hycom_nanflag' >> batch_files/file_${fileref}.ncl
              # Add time
               echo 'time = x1&time' >> batch_files/file_${fileref}.ncl
               # Write output file
               echo 'xxxx1!0 = "time"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time = time' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@standard_name = "time"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@long_name_name = "time"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@units = "days since 1850-01-01 00:00:00"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&time@calendar = "standard"'  >> batch_files/file_${fileref}.ncl
               echo 'xxxx1!1 = "lat"' >> batch_files/file_${fileref}.ncl
               # Add latitude
               echo 'xxxx1&lat = yoo' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@data_type = "float"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@units = "degrees_north"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@axis = "Y"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@standard_name = "latitude"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@valid_min = -90.0' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lat@valid_max =  90.0' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1!2 = "lon"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@data_type = "float"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@units = "degrees"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@axis = "X"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@standard_name = "longitude"' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@valid_min = -180.' >> batch_files/file_${fileref}.ncl
               echo 'xxxx1&lon@valid_max = 179.5' >> batch_files/file_${fileref}.ncl
               echo 'fout->water_v = xxxx1' >> batch_files/file_${fileref}.ncl

               rm -rf batch_files/file_${fileref}_submit.sh
               echo "#BSUB -J serialjob" > batch_files/file_${fileref}_submit.sh
               chmod u+x batch_files/file_${fileref}_submit.sh

               # Submit script
               echo "#BSUB -o logs/%J.out" >> batch_files/file_${fileref}_submit.sh
               echo "#BSUB -e logs/%J.err" >> batch_files/file_${fileref}_submit.sh
               echo "#BSUB -W 0:02" >> batch_files/file_${fileref}_submit.sh
               echo "#BSUB -q general" >> batch_files/file_${fileref}_submit.sh
               echo "#BSUB -n 1" >> batch_files/file_${fileref}_submit.sh
               echo "#" >> batch_files/file_${fileref}_submit.sh
               echo "ncl batch_files/file_${fileref}.ncl" >> batch_files/file_${fileref}_submit.sh

               # Check if out file exists
               if [[ ! -f ${odir}hycom_u_GLBu_${expnum}_${year}${month}${day}00_t0${hour}_sfc10m.nc ]]; then

                  echo "creating file ${odir}hycom_u_GLBu_${expnum}_${year}${month}${day}00_t0${hour}_sfc10m.nc"
                  bsub < batch_files/file_${fileref}_submit.sh
                  gettingafile=1
                  let fileref=$fileref+1
               else
                  echo "file ${odir}hycom_u_GLBu_${expnum}_${year}${month}${day}00_t0${hour}_sfc10m.nc exists"
                  if [[ ${gettingafile} -ne 1 ]];then
                     gettingafile=0
                  fi
               fi

            fi
         # End hour loop
         done
      # End day loop
      done
   # End month loop
   done
# End year loop
done
