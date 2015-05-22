function [flowu] = flowu(ucurve,vfir,vsec,cint,cslope,cage,wm,wf,nt,at,dmt,dft,y,t)
%wm,wf,nt,at as in the paper.
%dmt=1 if male works, dft=1 if female works.
%y outside income, t age.
if dmt==1&&dft==1
flowu=((wm+wf+y)^(1-ucurve))/(1-ucurve)+vfir*nt-vsec*nt^2-max(0,cint-cslope.*at-cage.*at.*t);
elseif dmt==1&&dft==0
flowu=((wm+y)^(1-ucurve))/(1-ucurve)+vfir*nt-vsec*nt^2;
elseif dmt==0&&dft==1
flowu=((wf+y)^(1-ucurve))/(1-ucurve)+vfir*nt-vsec*nt^2;
else
flowu=(y^(1-ucurve))/(1-ucurve)+vfir*nt-vsec*nt^2;
end

