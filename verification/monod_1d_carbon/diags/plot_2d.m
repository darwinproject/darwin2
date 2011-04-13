clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE SOME CHOICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtype=1;  %1=1d-k; 1.5=1.5d-j; 2=2d-jk
kn=1;       % for 1.5d,2d chose the level toplot
typlot=[1 2 3 7 8 9]; %1=t/par; 2=nutrients, 3=phy/zoo, 4=dom, 5=pom, 6=Rstar
                      % 7=carbon; 8=o2 and alk; 9=surface flux,ph, 10=Chl
dplot=2;  % 1 or 2d

geider=0;

rnp=16;
rsip=16;
rfep=1e-3;
rcp=120;
rop=170;

%%%%%% chose the run directory/time
rundir=['../run/'];
outputdir=['dar1d_20110407_0002/'];

timestr=['0000000000'];

npmax=9;
nzmax=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read the model output 
ncload([rundir,outputdir,'dic_tave.',timestr,'.t001.nc']);
ncload([rundir,outputdir,'ptr_tave.',timestr,'.t001.nc']);
ncload([rundir,outputdir,'dar_tave.',timestr,'.t001.nc']);
ncload([rundir,outputdir,'tave.',timestr,'.t001.nc']);
%%%%%%%%%%%%%%%%%%%

if (geider==0),
eval(['load ',rundir,'plankton_ini_char_nohead.dat']);
plankton=plankton_ini_char_nohead;
else
eval(['load ',rundir,'gplankton_ini_char_nohead.dat']);
plankton=gplankton_ini_char_nohead;
end


if (dtype==1);   yl=Z;  cay1=200; cay2=0;     end
if (dtype==1.5); xl=Y;  cax1=1; cax2=65;       end


%%%%%%%%%%%%%%%%% MAIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%
if (dplot==2),
 for itype=typlot,
%%%%%%%%%%
  ncon=5;
  if (itype==1)
% plot T/PAR
   spmax=2;
   tmp1=Ttave; min1=4; max1=12; tstr1='Theta (C)';  np1=1;
   tmp2=dar_PAR_ave; min2=0; max2=220;
                                 tstr2='PAR (\mu Ein/m2/s)'; np2=1;
   ncon=10;
  end % if itype 1
  if (itype==2)
% plot nutrients
   spmax=4;
   tmp1=PO4; min1=0; max1=2; tstr1='PO4 (\muM)';  np1=1;
   tmp2=NO3; min2=0; max2=20; tstr2='NO3 (\muM)';  np2=1;
   tmp3=FeT*1e3; min3=0; max3=1.4; tstr3='Fe (nM)';  np3=1;
   tmp4=SiO2; min4=0; max4=10; tstr4='Si (\muM)';  np4=1;
    ncon=5;
  end % if itype 2
  if (itype==3)
% plot plankton
   spmax=3;
   min1=0; max1=.1; tstr1='Phyto biomass (\muM P)';  np1=npmax;
   min3=0; max3=1e-2; tstr3='Zoo biomass (\muM P)';  np3=nzmax;
   min2=0; max2=1; tstr2='Chl (mg/m3)'; np2=npmax;
  end % if itype 3
  if (itype==4)
% plot DOM
   spmax=3;
   tmp1=DOP; min1=0; max1=0.1; tstr1='DOP (\muM)';  np1=1;
   tmp2=DON; min2=0; max2=3; tstr2='DON (\muM)';  np2=1;
   tmp3=DOFe*1e3; min3=0; max3=0.4; tstr3='DOFe (nM)';  np3=1;
  end % if itype 4
  if (itype==5)
% plot POM
   spmax=4;
   tmp1=POP; min1=0; max1=0.1; tstr1='POP (\muM)';  np1=1;
   tmp2=PON; min2=0; max2=.2; tstr2='PON (\muM)';  np2=1;
   tmp3=POFe*1e3; min3=0; max3=0.1; tstr3='POFe (nM)';  np3=1;
   tmp4=POSi; min4=0; max4=2; tstr4='POSi (\muM)';  np4=1;
  end % if itype 3
%%%%%%%%%
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
   tmp1=DIC; min1=1800; max1=2200; tstr1='DIC (\muM)';  np1=1;
   tmp2=DOC; min2=0; max2=50; tstr2='DOC (\muM)';  np2=1;
   tmp3=POC; min3=0; max3=5; tstr3='POC (nM)';  np3=1;
   tmp4=PIC; min4=0; max4=2; tstr4='PIC (\muM)';  np4=1;
  end % if itype 3
  if (itype==8)
% plot o2 alk
   spmax=2;
   tmp1=ALK; min1=2250; max1=2350; tstr1='ALK (\muM eq)';  np1=1;
   tmp2=O2; min2=0; max2=500; tstr2='O2 (\muM)';  np2=1;
  end % if itype 6
  if (itype==9)
% plot fluxes
   spmax=4;
   tmp1=dic_SUR_ave; min1=-4; max1=4; tstr1='CO2 Flux (mmol/m2/s)';  np1=1;
   tmp2=dic_SURO_ave; min2=-4; max2=4; tstr2='O2 Flux (mmol/,2/s)';  np2=1;
   tmp3=dic_pH_ave; min1=7.8; max1=8.2; tstr3='pH';  np1=1;
   tmp4=dic_pCO2_ave; min2=260e-6; max2=300e-6; tstr4='pCO2';  np2=1;
  end % if itype 9
% plot Chl
  if (itype==10)
% plot chl
   spmax=1;
   tmp1=dar_Chl_ave; min1=0; max1=1; tstr1='Chl (mg/m3)';  np1=1;
  end % if itype 10




%% Plot time/depth output

    figure(itype)
    clf

    if (itype~=3 & itype~=9) % loop through figures
        for isp=1:spmax % loop through subplots
            subplot(spmax,1,isp)

            eval(['plottmp=tmp' num2str(isp) ';']);
%           contourf(T./31104000,[0;-Z(1:end-1)],plottmp',ncon),axis ij
            pcolor(T./31104000,[0;-Z(1:end-1)],plottmp'),axis ij
            shading flat, 
            ylim([0 cay1])

            eval(['title(tstr' num2str(isp) ')'])
            eval(['caxis([min' num2str(isp) ',max' num2str(isp) '])'])
            colorbar
            colormap(flipud(hot))
        end
    end
   if itype==3 % loop through figures
        for isp=1:spmax % loop through subplots
            eval(['npn=np',num2str(isp),';']);


            if isp==1; ptype='Phy'; tag=''; end
            if isp==3; ptype='ZOO'; tag='P'; end
%           if isp==2; ptype='Chl'; tag=''; end
            for ip=1:npn % sum all phyto/zoo-plankton types
                if ip==1;
                    eval(['planktmp=zeros(size(' ptype '1' tag '));']);
                end
                eval(['planktmp=planktmp+' ptype num2str(ip) tag ';']);
            end
            subplot(spmax,1,isp)

            %contourf(T./31104000,[0;-Z(1:end-1)],planktmp',ncon),axis ij
            pcolor(T./31104000,[0;-Z(1:end-1)],planktmp'),axis ij
            shading flat, 
            ylim([0 cay1])

            eval(['title(tstr' num2str(isp) ')'])
            eval(['caxis([min' num2str(isp) ',max' num2str(isp) '])'])
            colorbar
            colormap(flipud(hot))
        end
    end
    if (itype==9) % loop through figures
        for isp=1:spmax % loop through subplots
            subplot(spmax,1,isp)

            eval(['plottmp=tmp' num2str(isp) ';']);
%           contourf(T./31104000,[0;-Z(1:end-1)],plottmp',ncon),axis ij
            plot(T./31104000,plottmp','k-'),

            eval(['title(tstr' num2str(isp) ')'])
%           eval(['caxis([min' num2str(isp) ',max' num2str(isp) '])'])
        end
    end

end

end % if dplot
