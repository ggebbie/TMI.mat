function [EDIC,ETALK] = get_E_CO3(E,cDIC,cTALK,pressure)

% My functional form of CO3 vs (TALK-DIC) didn't work for the
% sensitivities.
% Instead try an impulse response method.
[iii,jjj,sss]=find(E);
Nfield = length(cDIC);
ipert = unique(jjj);
sensDIC  = zeros(length(ipert),1);
sensTALK = zeros(length(ipert),1);
for ip = 1:length(ipert)
  delta = 1e-3;
  tmp0 =  CO3_from_DIC_ALK(cDIC(ipert(ip)),cTALK(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
  tmp1 =  CO3_from_DIC_ALK(cDIC(ipert(ip))+delta,cTALK(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
  tmp2 =  CO3_from_DIC_ALK(cDIC(ipert(ip)),cTALK(ipert(ip))+delta,35,10,pressure(ipert(ip)),50,1,[]); 
  sensDIC(ip) = (tmp1-tmp0)./delta;
  sensTALK(ip) = (tmp2-tmp0)./delta;
  
  %senstest(ip) = -bb-2.*cc.*cCO3(ipert(ip));
end
Etmp = sparse(ipert,ipert,sensDIC,Nfield,Nfield);
EDIC = E*Etmp;

Etmp = sparse(ipert,ipert,sensTALK,Nfield,Nfield);
ETALK = E*Etmp;

