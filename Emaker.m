function [E] = Emaker(obslon,obslat,obsdep,gridlon,gridlat,kgrid, ...
                      griddep,radius)
%function [E] = Emaker(obslon,obslat,obsdep,gridlon,gridlat,kgrid,griddep,radius)
%
% Example:
% ENd =  = Emaker(Ndlon,Ndlat,Nddep,grd.LON(grd.it),grd.LAT(grd.jt),grd.kt,grd.DEPTH,3.9)

r2 = radius.^2;
Nobs = length(obslon);
Ngrid = length(gridlon);
E = sparse(Nobs,Ngrid);
maxz = max(griddep);
Nz = length(griddep);
eps = 1e-4;

for nd = 1:Nobs
  nd;

  %% get two closest depths.
  level = find(griddep > obsdep(nd));  
  if isempty(level)
    level = Nz;
  end
  level = level(1);

  zweight(1) = 1./(abs(griddep(level)-obsdep(nd))+eps);
  zweight(2) = 1./(abs(griddep(level-1)-obsdep(nd))+eps);

  zweight = zweight./sum(zweight);

  %% two levels. do them separately.
  good2 = [];
  goodloc2 = []; 
  for nz = 1:2
    zlev = find(kgrid == level);

    dist = (gridlon(zlev)-obslon(nd)).^2 + (gridlat(zlev)-obslat(nd)).^2 ...
           + eps;
    if obslon(nd) >= 360-radius; 
      dist = min(dist,(gridlon(zlev)-obslon(nd)-360).^2  ...
                    + (gridlat(zlev)-obslat(nd)).^2) + eps;
    elseif obslon(nd) <= radius
      dist = min(dist,(gridlon(zlev)-obslon(nd)+360).^2  ... 
            + (gridlat(zlev)-obslat(nd)).^2) + eps;
    end
    %[sd,loc]=sort(dist);
    good = dist(dist<r2);
    goodloc = dist<r2;
    if ~isempty(good)
        %[good,goodloc] = min(dist);
        %end
      good2 = [good2; good./zweight(nz)];
      goodloc2 = [goodloc2; zlev(goodloc)];
    end
    level = level -1;
  end
  if~isempty(goodloc2)
    E(nd,goodloc2) = (1./good2)./sum(1./good2);
  end
end



