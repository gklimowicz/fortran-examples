#!/usr/bin/ksh
USAGE="$0  arg1 ..."
# SCRIPT: plot_hyc_msf.ksh
#
# PURPOSE: ploting script for meridional stream functions
#          prepared avg.f
# INPUT    Define number of files (1,2,3,4) 
#          and abs path input files
# OUTPUT   post script file 
#
# set -n   # Uncomment to check your syntax, without execution.
#          # NOTE: Do not forget to put the comment back in or
#          #       the shell script will not execute!
# set -x   # Uncomment to debug this shell script (Korn shell only)
#          
##########################################################
########### DEFINE FILES AND VARIABLES HERE ##############
##########################################################

#INPUT FILES:
num_files=1
pathname_1=/discover/nobackup/ntausnev/RUNS_ME/EhV2a/z_00HYC_diag/avg_ov_EhV2a_2100-2109.txt
pathname_2=/discover/nobackup/ntausnev/RUNS_ME/EhV2bh/z_00HYC_diag/avg_ov_EhV2bh_2100-2109.txt

pathname_3=/discover/nobackup/ntausnev/RUNS_ME/EhCld0/z_00HYC/avg_ov_EhCld0_239.txt
pathname_4=/discover/nobackup/ntausnev/RUNS_ME/EhMay2/z_00HYC/avg_ov_EhMay2_239.txt

#### You do not need to change the code below ###

filename_1=${pathname_1##*/}
filename_1_no_suffix=${filename_1%.txt}
filename_2=${pathname_2##*/}
filename_2_no_suffix=${filename_2%.txt}
filename_3=${pathname_3##*/}
filename_3_no_suffix=${filename_3%.txt}
filename_4=${pathname_4##*/}
filename_4_no_suffix=${filename_4%.txt}


#OUTPUT FILE(S):
pathname_eps=${filename_1_no_suffix}.eps

##########################################################
############### DEFINE FUNCTIONS HERE ####################
##########################################################

function help_use
{
#
# Display help message and quit
#
cat << ENDOFTEXT

Dear $USER , the usage of the script is as follows:
     $0 

ENDOFTEXT
exit 1
}

# Shun Sun matlab script:
cat > matlab_plot_msf.m << END_MATLAB_SCRIPT

% matlab  -nojvm -nosplash -nodisplay < matlab_plot_msf.m  >/dev/null 2>&1

clear; clf; 
clr=1;             % 1 for color; 0 for BW
add_bar=0;         % 1 for adding colorbar;

ii=387; kk=27; kdm=[1:kk];
load latlon387.txt; lat=latlon387(:,2);
%sig=[28.89 30.07 31.11 32.02 32.81 33.49 34.07 34.56 34.97 35.31 35.59 35.82 36.01 36.17 36.31 36.44 36.56 36.67 36.77 36.86 36.94 37.01 37.07 37.12 37.16 37.20];
% sig1g
sig=[24.35 26.07 27.50 28.67 29.61 30.35 30.92 31.35 31.67 31.90 32.06 32.17 32.25 32.31 32.36 32.40 32.43 32.46 32.49 32.52 32.54 32.56 32.58 32.60 32.62 32.64];
          ]

imin=[138 203 140 138];imax=[293 293 299 360]; %for AT,IN,PA,GL
ia=imin;ib=imax;

vv=[-50:2:50];
v1=[-50:4:-2]; v2=[2:4:34]; vv1=[-50:8:-2];vv2=[2:8:42]; xmin=-75; xmax=70;
v1=[-50:2:-2]; v2=[2:2:34]; vv1=[-50:4:-2];vv2=[2:4:42]; xmin=-75; xmax=70;
vspace=0.05; hspace=.005;

%nnfile=1;  % number of input data      % TNL
nnfile=$num_files;  % number of input data      % TNL

j1=.99;
for nfile=1:nnfile
w=.50;hgt=.34;i1=.06;i2=i1+hspace+w; 
if (nfile==nnfile) hgt=.50; vspace=0.07; end
j1=j1-hgt; hgt=hgt-vspace; dx=.00145;

clear ov4;
if nfile==1
 % load /discover/nobackup/ntausnev/RUNS_ME/EhCld0/z_00HYC/avg_ov_EhCld0_239.txt; ov4=avg_ov_EhCld0_239; % Original
 % load ${pathname_1}; ov4=${filename_1_no_suffix}; % TNL
 ov4=load('${pathname_1}'); 
elseif nfile==2
 ov4=load('${pathname_2}'); 
elseif nfile==3
 ov4=load('${pathname_3}'); 
elseif nfile==4
 ov4=load('${pathname_4}'); 
end

for n=1:4
w1=dx*(ib(n)-ia(n));
if n==1
  axes('position',[i1,j1,w1,hgt]);
elseif n==2
  axes('position',[i1,j1,w1,hgt]);
elseif n==3
  axes('position',[i1,j1,w1,hgt]);
elseif n==4
  axes('position',[i1,j1,w1,hgt]);
if (clr==1) & (n==4) & (nfile==nnfile) & (add_bar==1)
  axes('position',[i2,j1-0.0475,w,hgt+0.0475]);     
end
end
i1=i1+hspace+w1;
ov1=ov4(1+(n-1)*(kk-1):n*(kk-1),:);
ov(1,1:ii)=0; ov(2:kk,1:ii)=ov1;

% --- mask out unwanted domain
for i=1:imin(n)-1
  ov(:,i)=nan;
end
for i=imax(n)+1:ii
  ov(:,i)=nan;
end
if clr==1
  [c,h]=contourf(lat',kdm',ov,vv); hold on
  caxis([-40 40]);
  for icnt = 1: length(h)
  set( h(icnt), 'EdgeColor', 'none' )
  end
end
[c,h]=contour (lat',kdm',ov,v1,'k--'); hold on
if h>0 clabel(c,h,vv1,'LabelSpacing',100,'FontSize',9)  % 72 is one inch
%if n==2 clabel(c,h,v1,'LabelSpacing',72,'FontSize',9)  % 72 is one inch
%else clabel(c,h,vv1,'LabelSpacing',72,'FontSize',9); end
end
[c,h]=contour (lat',kdm',ov,v2,'k'); hold on
if h>0 clabel(c,h,vv2,'LabelSpacing',100,'FontSize',9)  % 72 is one inch
%if n==2 clabel(c,h,v2,'LabelSpacing',72,'FontSize',9)  % 72 is one inch
%else clabel(c,h,vv2,'LabelSpacing',72,'FontSize',9); end
end
set(gca,'XTick',[-60:30:60]);
%set(gca,'YTicklabel',sig)
set(gca,'YTick',[1:1:kk]);
axis ij;
if n==1 
   if nfile==1   title ('Atlantic','Fontsize',13);
%elseif nfile==2 title ('Atlantic','Fontsize',13); 
%elseif nfile==3 title 'a20t3 yr230 Atlantic'; 
%elseif nfile==4 title 'a20t6 yr20 Atlantic'; 
  end
%set(gca,'XTicklabel',['-60';'-40';'-20';'  0';' 20';' 40';'   '],'FontSize',8);
%text(-67,kdm(3),'Atlantic','color','k','Fontsize',9)
text(lat(ib(n))-20,kk/2+1,'\sigma', 'color','k','rotation',90,'fontsize',9);
  for i=1.5:kk-.5
  text(lat(ib(n))-12,i,num2str(sig(i-.5),'%5.2f'),'color','k','Fontsize',7)
  end
elseif n==2 & nfile==1
     title ('Indian','Fontsize',13); 
%text(-67,kdm(3),'Indian','color','k','Fontsize',9)
elseif n==3 & nfile==1
     title ('Pacific','Fontsize',13); 
%text(-67,kdm(3),'Pacific','color','k','Fontsize',9)
     if (nfile==nnfile) xlabel 'Latitude';end
elseif n==4 & nfile==1
     title ('Global','Fontsize',13); 
%text(-67,kdm(3),'Global','color','k','Fontsize',9)
end
set(gca,'YTicklabel',['';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';''],'FontSize',8);
axis ([lat(ib(n)) lat(ia(n)) 1 kk])

if (clr==1) & (n==4) & (nfile==nnfile) & (add_bar==1)
  h=colorbar('horiz');
  set(h,'Position',[0.39 0.05 0.3 0.02]);             %line 77
end                                                 %line 78

end  % for n=1:4

end   % for nfile=1:nnfile

orient landscape
% print -depsc2 ov_1deg.eps      % TNL original
print -depsc2 ${pathname_eps}
%print -djpeg99 ov_1deg.jpg      % TNL comment

%%%%%%%%%%%% END_MATLAB_SCRIPT   %%%%%%%%%%%%%%%%%%%
END_MATLAB_SCRIPT

##########################################################
################ BEGINNING OF MAIN #######################
##########################################################
 
print "##################################################################################"
print "Script: $0  starts: ==========> "                 
if (( $# != 0 )) 
then
    help_use
fi

matlab  -nojvm -nosplash -nodisplay < matlab_plot_msf.m  >/dev/null 2>&1

cat << ENDOFTEXT
Script: $0 ended
output file is: ${pathname_eps}
+++++ Normal end of execution +++++
##################################################################################
ENDOFTEXT
# End of script

