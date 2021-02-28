
%%
filenc = 'tracer_output.nc';
ncid = netcdf.create(filenc,'NC_WRITE');

dimid_time = netcdf.defDim(ncid,'TIME',length(years))

dimid_depth = netcdf.defDim(ncid,'DEPTH',max(kt))
dimid_lat = netcdf.defDim(ncid,'LAT',max(jt))
dimid_lon = netcdf.defDim(ncid,'LON',max(it))
%
varid_lon = netcdf.defVar(ncid,'longitude','double',dimid_lon);
netcdf.putAtt(ncid,varid_lon,'long_name','Longitude')
netcdf.putAtt(ncid,varid_lon,'units','degrees_east')
% 
varid_lat = netcdf.defVar(ncid,'latitude','double',dimid_lat);
netcdf.putAtt(ncid,varid_lat,'long_name','Latitude')
netcdf.putAtt(ncid,varid_lat,'units','degrees_north')
%
varid_depth = netcdf.defVar(ncid,'depth','double',dimid_depth);
netcdf.putAtt(ncid,varid_depth,'long_name','Depth')
netcdf.putAtt(ncid,varid_depth,'units','meters')
%
varid_tracer = netcdf.defVar(ncid,'C','double', ...
                       [dimid_time dimid_depth dimid_lat dimid_lon]);
netcdf.putAtt(ncid,varid_tracer,'long_name','passive tracer/temperature/d18O')
netcdf.putAtt(ncid,varid_tracer,'units','[C-stuff/kg]')
%
%
varid_notes = netcdf.defVar(ncid,'notes','double',0);
netcdf.putAtt(ncid,varid_notes,'author','Jake Gebbie (jgebbie@whoi.edu)')
netcdf.putAtt(ncid,varid_notes,'version','Gebbie Huybers 2012 transient driver')
netcdf.putAtt(ncid,varid_notes,'output created','22 Sep 2018')
netcdf.putAtt(ncid,varid_notes,'file created','22 Sep 2018')
%
netcdf.endDef(ncid)
% % % Agregar datos de coordenadas
netcdf.putVar(ncid,varid_lon,LON);
netcdf.putVar(ncid,varid_lat,LAT);
netcdf.putVar(ncid,varid_depth,DEPTH);
netcdf.putVar(ncid,varid_tracer,Cfield);
netcdf.close(ncid) 



