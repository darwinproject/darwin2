clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE SOME CHOICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtype=1;  %1=1d-k; 1.5=1.5d-j; 2=2d-jk
kn=1;       % for 1.5d,2d chose the level toplot
typlot=[3 ]; %1=t/par; 2=nutrients, 3=phy/zoo, 4=dom, 5=pom, 6=Rstar
                      % 7=carbon % 8=o2 and alk
rnp=16;
rsip=16;
rfep=1e-3;
rcp=120;
rop=170;

%%%%%% chose the run directory/time
rundir=['../run/'];
outputdir=['dar1d_20090609_0001/'];

timestr=['0000000000'];

npmax=2;
nzmax=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read the model output 
ncload([rundir,outputdir,'ptr_tave.',timestr,'.t001.nc']);
ncload([rundir,outputdir,'dar_tave.',timestr,'.t001.nc']);
ncload([rundir,outputdir,'tave.',timestr,'.t001.nc']);
%%%%%%%%%%%%%%%%%%%

eval(['load ',rundir,'plankton_ini_char_nohead.dat']);
plankton=plankton_ini_char_nohead;


if (dtype==1);   yl=Z;  cay1=-200; cay2=0;     end
if (dtype==1.5); xl=Y;  cax1=1; cax2=65;       end


%%%%%%%%%%%%%%%%% MAIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%
for it=1:length(T)
 for itype=typlot,
%%%%%%%%%%
  if (itype==1)
% plot T/PAR
   spmax=2;
   tmp1=Ttave; min1=4; max1=12; tstr1='Theta (C)';  np1=1;
   tmp2=dar_PAR_ave; min2=0; max2=220; 
                                 tstr2='PAR (\mu Ein/m2/s)'; np2=1;
  end % if itype 1
  if (itype==2)
% plot nutrients
   spmax=4;
   tmp1=PO4; min1=0; max1=2; tstr1='PO4 (\muM)';  np1=1;
   tmp2=NO3; min2=0; max2=20; tstr2='NO3 (\muM)';  np2=1;
   tmp3=FeT*1e3; min3=0; max3=1.4; tstr3='Fe (nM)';  np3=1;
   tmp4=SiO2; min4=0; max4=20; tstr4='Si (\muM)';  np4=1;
  end % if itype 2
  if (itype==3)
% plot plankton
   spmax=2;
   min1=0; max1=1e-1; tstr1='Phyto biomass (\muM P)';  np1=npmax;
   min2=0; max2=1e-2; tstr2='Zoo biomass (\muM P)';  np2=nzmax;
  end % if itype 3
  if (itype==4)
% plot DOM
   spmax=3;
   tmp1=DOP; min1=0; max1=0.1; tstr1='DOP (\muM)';  np1=1;
   tmp2=DON; min2=0; max2=3; tstr2='DON (\muM)';  np2=1;
   tmp3=DOFe*1e3; min3=0; max3=0.4; tstr3='DOFe (nM)';  np3=1;
  end % if itype 4
  if (itype==5)
% plot plankton
   spmax=4;
   tmp1=POP; min1=0; max1=0.1; tstr1='POP (\muM)';  np1=1;
   tmp2=PON; min2=0; max2=.2; tstr2='PON (\muM)';  np2=1;
   tmp3=POFe*1e3; min3=0; max3=0.1; tstr3='POFe (nM)';  np3=1;
   tmp4=POSi; min4=0; max4=2; tstr4='POSi (\muM)';  np4=1;
  end % if itype 3
%%%%%%%%%
  if (itype==6),
   spmax=2;
   min1=0; max1=.05; tstr1='Rstar (\muM P)';  np1=npmax;
   tmp2=PO4; min2=0; max2=0.05; tstr2='PO4 (\muM)';  np2=1;
  end % if itype 6
%%%
  if (itype==7)
% plot carbon
   spmax=4;
   tmp1=DIC; min1=1000; max1=2500; tstr1='DIC (\muM)';  np1=1;
   tmp2=DOC; min2=0; max2=5; tstr2='DOC (\muM)';  np2=1;
   tmp3=POC; min3=0; max3=2; tstr3='POC (nM)';  np3=1;
   tmp4=PIC; min4=0; max4=2; tstr4='PIC (\muM)';  np4=1;
  end % if itype 3
  if (itype==8)
% plot o2 alk
   spmax=2;
   tmp1=ALK; min1=1000; max1=2800; tstr1='ALK (\muM eq)';  np1=1;
   tmp2=O2; min2=0; max2=400; tstr2='O2 (\muM)';  np2=1;
  end % if itype 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(itype), if (it==1), hold off, subplot(1,1,1); end
  if (dtype==1),  spl1=1; spl2=spmax; end
  if (dtype==1.5), spl1=spmax; spl2=1; end
  for isp=1:spmax,
   if (dtype==1),  
         eval(['cax1=min',num2str(isp),';']);
         eval(['cax2=max',num2str(isp),';']);
   end
   if (dtype==1.5),
         eval(['cay1=min',num2str(isp),';']);
         eval(['cay2=max',num2str(isp),';']);
   end
   eval(['npn=np',num2str(isp),';']);
   eval(['tstr=tstr',num2str(isp),';']);
   subplot(spl1,spl2,isp); hold off
%%%%%%%%%%%%%%%%%%%%%%
   ip1=0;
   for ip=1:npn,
    if (itype==3 | (itype==6 & isp==1) ), 
       ip1=ip1+1; if (ip1==8), ip1=1; end
       clear tmp,
       if (itype==3),
        if (isp==1), 
          if (ip<10), eval(['tmp=Phy',num2str(ip),';']); end
          if (ip>9), eval(['tmp=Phy',num2str(ip),';']); end
        end
        if (isp==2), eval(['tmp=ZOO',num2str(ip),'P;']); end
       end % if
       if (itype==6),
        if (isp==1), eval(['tmp=dar_Rstar_Phy',num2str(ip),';']); end
       end % if
       if (isp==1)
         if (plankton(ip,3)==1), 
           if (plankton(ip,1)==1), 
             cstr=['r-']; 
           end
           if (plankton(ip,1)==2),
             cstr=['y-'];
           end
           if (plankton(ip,1)==0),
             cstr=['k-'];
           end % if
         else
           if (plankton(ip,19)==3),
             cstr=['b-'];
           else
             cstr=['g-']; 
           end % if
         end % if
       end % if
       if (isp==2), 
        if (ip==1), cstr=['k-']; end
        if (ip==2), cstr=['r--']; end
        if (ip==3), cstr=['g:']; end
       end % if
%      if (ip1==1), cstr=['k-']; end
%      if (ip1==2), cstr=['r-']; end
%      if (ip1==3), cstr=['b-']; end
%      if (ip1==4), cstr=['g-']; end
%      if (ip1==5), cstr=['y-']; end
%      if (ip1==6), cstr=['c-']; end
%      if (ip1==7), cstr=['m-']; end
     else
        clear tmp, eval(['tmp=tmp',num2str(isp),';']);
        cstr='b-';
     end % if not itype 3
%%%
     if (dtype==1), xl=squeeze(tmp(it,:)); end
     if (dtype==1.5), yl=squeeze(tmp(it,kn,:)); end
     if (itype==6),
        clear fi, fi=find(yl<0); yl(fi)=NaN;
        clear fi, fi=find(isinf(yl)==1); yl(fi)=NaN;
     end % if
     plot(xl,yl,cstr); hold on, axis([cax1 cax2 cay1 cay2]);
   end % for ipn
%%%%%%%%%%%%%%%%%%%%%%
   if (dtype==1), xlabel(tstr); end
   if (dtype==1.5), ylabel(tstr); end 
   if (isp==1), title(['Time=',num2str(it)]); end
  end % for isp
%%%%%%%%%%%%%%%%
 end % for itype
 pause(2)
end % for it
