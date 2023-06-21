%script intracellular distribution of molecule.m
%Author: KH.
%Varshney Lab, OMRF, OKC, OK.
clear all %Make sure not other scripts are running before execution of this script.

%load packages for analysis and define contatins function which is absent in GNU Octave 
pkg load signal
contains=@(str, pattern)~cellfun('isempty',strfind(str, pattern))
%instruction box
uiwait(msgbox("Please make sure to have the cell contour marker intensity profile and the intracellular molecule intensity profile included as separate .csv files in you current directory. \nFor each intensity profile files, make sure each column represent the intensity profile readout of individual sample. \nThe following prompt will ask for the name of cell contour marker and intracellular molecule which you should make sure is included as part of the file name. \ni.e. cell contour marker is HCS1 and the file name is wt hcs1.csv, enter in prompt 'hcs1'. \nDo note that this analysis is designed for only 2-dimension and for mapping only one molecule."));
%------------------get response from input for cell contour-----------------------
response=inputdlg({'Cell contour','intracellular molecule'});
contour=response(1)
probe=response(2)
%------------------get resolution of images---------------------------------------
imageresolution=inputdlg({'Resolution of the image','Unit of resolution'});
resolution=imageresolution(1)
imageunit=imageresolution(2)
%---------------------------------------------------------------------------------
working_directory=pwd;
files=dir('*.csv');
filenames={files.name};
contourfilematch=contains(filenames, contour)
probefilematch=contains(filenames, probe)
if any(contourfilematch);
  contournum=find(contourfilematch==1)
  c=filenames(contournum)
  contourprofile=csvread(char(c))
  else
errordlg('Files not found in directory');
end
if any(probefilematch);
  probenum=find(probefilematch==1)
  p=filenames(probenum)
  probeprofile=csvread(char(p))
else
errordlg('Files not found in directory');
end
contourprofilesmoothed=abs(sgolayfilt(contourprofile));
[nrow ncol]=size(contourprofilesmoothed);
probeintensity=probeprofile;
contourmap=zeros(nrow, ncol);
probemap=zeros(nrow,ncol);
for n=1:ncol
contourpeak=strcat('contour_sample',num2str(n),'.mat');
[conval conloc conextra]= findpeaks(contourprofilesmoothed(1:end,n));
save(contourpeak, "conloc", "conval", "conextra");
probepeak=strcat('probe_sample',num2str(n),'.mat');
[proval proloc proextra]=findpeaks(probeprofile(1:end,n));
save(probepeak, "proloc", "proval", "proextra");
end
clear -v conloc conval conextra proval proloc proextra n
%following the identification of the contour peaks, find the leftmost tallest peak and its x and y value for x-shift superimpose
matfiles=dir('*.mat');
matfilenames={matfiles.name};
contourprofilenormalized=contourprofilesmoothed;
for n=1:ncol
  samplesearch=strcat('sample',num2str(n),'.mat');
  matsearchresult=contains(matfilenames, samplesearch);
  toopen=matfilenames(find(matsearchresult==1));
  % This creates a list of .mat files to open after creating the contour and probe mat files. Since the files are established in with contour and probe up front and as such are sorted alphabetically, additional sorting rule is not needed for this script.
  %contour is always toopen{1} and provides conval and conloc variables and probe is always toopen{2} and provides proval and proloc variables.
  load(toopen{1});
  load(toopen{2});
maxcontourintensity=max(conval);
contourprocessing=(contourprofilesmoothed(1:end,n)/maxcontourintensity);
maxprobeintensity=max(proval);
probeprocessing=(probeintensity(1:end,n)/maxprobeintensity);
  for idx=n
  contourprofilenormalized(1:end,idx)=contourprocessing;
  %use normalized profile to find the peaks again for cell base alignment.
  [val loc]=findpeaks(contourprofilenormalized(1:end,idx), "MinPeakDistance", 50);
  %Find the position within array the max value of the peak.
  [maxvalrow maxvalcol]=max(val)
  %Examine whether the max peak is the left-most major max peak (base start from 0 so defined base of cell as left-most major max).
  default_contourbase_value=loc(maxvalcol)
  sortval=sort(val)
    %sort value from lowest to highest, so we can determine the second highest peak.
  searchresult= find(val==(sortval(end-1)))
    %determine the position of the second highest peak in respect to the val rows. 
  secondpeakloc=loc(searchresult)
    %the loc position of the second highest peak.
    
  if secondpeakloc<default_contourbase_value
    %if the second tallest peak is the left of the tallest peak, then the second tallest peak is the base and should be used as the reference point.
    shiftx=secondpeakloc
  else
    shiftx=default_contourbase_value;
  end
arrayrownum=length(contourprofilenormalized(shiftx:end,idx));
contourmaptoadd=(contourprofilenormalized(shiftx:end,idx));
probemaptoadd=(probeprocessing(shiftx:end));
end
contourmap(1:arrayrownum,n)=contourmaptoadd;
probemap(1:arrayrownum,n)=probemaptoadd;
end
%Statistical distribution
statscontour=statistics(contourmap,2);
averagecontour=statscontour(:,5);
errcontour=statscontour(:,6);
csvwrite('statistical analysis of contour.txt',statscontour);
statsprobe=statistics(probemap,2);
averageprobe=statsprobe(:,5);
errprobe=statsprobe(:,6);
csvwrite('statistical analysis of probe.txt',statsprobe);
%----------------map for probe intensity > normalized 51%-----------------
probelength=probemap>0.51
%----------------Make figures---------------------------------------------
xlabelcontent=strcat('Position(',imageunit,')');
xmax=(str2num(resolution{1}))*nrow;
xint=(str2num(resolution{1}));
x=0+(0:nrow-1)*xint;
figure(1);
plot(x, averagecontour,'m');
hold on;
plot(x, averageprobe,'g');
xlabel(xlabelcontent);
ylabel('Fluorescent intensity proportion');
ylim([0 1.2]);
grid on;
set(gca,'Color','k');
set(gcf,'InvertHardCopy','off');
File_Name1=inputdlg({'Please name your figure to save, the file will be saved as .png'});
saveas(figure(1),File_Name1{1},'png');
%Load CSV data file of cell contour intensity as measured in ImageJ.