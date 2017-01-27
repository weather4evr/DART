#!/bin/csh

set BASE = /glade/p/work/thoar/roms/test/WC13/OldData/ 

foreach FILE ( obs_wc13_merged_2013+07d_pp1_depthinmeters_dt5400_nobio.nc \
               obs_wc13_merged_2013+07d_pp1_depthinmeters_dt5400_nobio.nc \
               obs_wc13_merged_2013+07d_pp1_depthinmeters_dt5400_nobio.nc )

cp -v ${BASE}/${FILE} .

ncatted -a    units,survey_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
        -a calendar,survey_time,c,c,'gregorian' \
        -a    units,obs_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
        -a calendar,obs_time,c,c,'gregorian' \
        ${FILE}

end