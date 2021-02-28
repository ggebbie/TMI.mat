function [CO3] = CO3_from_DIC_ALK(DIC,ALK,SAL,TempC,P,Si,PO4,options)
% function [CO3] = CO3_from_DIC_ALK(DIC,ALK,SAL,TempC,P,Si,PO4,options)
%
% CO3_from_DIC_ALK        Carbonate ion concentration from
%                         dissolved inorganic carbon (DIC) and
%                         total alkalinity (ALK)
%==========================================================================
%
% USAGE:
% [CO3] = CO3_from_DIC_ALK(DIC,ALK,SAL,TempC,P,Si,PO4,{options}) 
%
% DESCRIPTION:
%  Calculates CO3(2-) carbonate ion concentration of seawater from dissolved inorganic carbon, total
%  alkalinity, salinity, temperature, pressure, silicate, and phosphate
%
% INPUT:
%  DIC =  Dissolved Inorganic Carbon                              [ umol/kg ]
%  ALK =  Total Alkalinity                                        [ umol/kg ]
%  SAL =  Practical Salinity                                      [ PSS-1978] 
%  TempC = In-situ temperature                                    [  deg C  ]
%  P     = Pressure                                               [  dbar   ]
%  Si    = silicate                                               [ umol/kg ]
%  PO4   = Phosphate                                              [ umol/kg ]
%  OPTIONS (optional) : pHscale = 1 = total scale
%           K1K2    = 'mehrbach-dickson-millero'
%           KSO4    = 'dickson-lee'
%
% OUTPUT:
%  CO3  =  Carbonate ion concentration                            [ umol/kg ]
%
% AUTHOR:
%  G. Jake Gebbie, WHOI
%  
% Adapted from CO2SYS:
%    Lewis, E., and D. W. R. Wallace. 1998. Program Developed for
%    CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information
%    Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy,
%    Oak Ridge, Tennessee. 
%    http://cdiac.ornl.gov/oceans/co2rprt.html
%  
% VERSION NUMBER: 0.1 (18 Nov 2015)
%
% REFERENCES:
%
%==========================================================================
%varargin
if isempty(options)
  % default options.
  options.pHscale = 'total'
  options.K1K2    = 'mehrbach-dickson-millero'
  options.KSO4    = 'dickson-lee'
end


TA=ALK./1e6; % Convert from micromol/kg to mol/kg
TC=DIC./1e6; % Convert from microatm. to atm.
TP=PO4./1e6;
TSi = Si./1e6;

[K0,fH,FugFac,VPFac,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,TB,TF,TS] = Constants(SAL,TempC,P,options);
%[K0,fH,FugFac,VPFac,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,TB,TF,TS] = Constants_orig(SAL,TempC,P);

pH   = pH_from_DIC_ALK(TC,TA,K1,K2,KW,KP1,KP2,KP3,TP,TSi,KSi,TB,KB,TS,KS,TF,KF);
CO3  = CO3_from_pH_DIC(pH,DIC,K1,K2);

