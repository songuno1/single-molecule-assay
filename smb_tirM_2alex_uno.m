function smb_tirM_2alex_uno
% Single Molecule Biophysics Lab. in Seoul National University 

clear all;  close all;  disp('           ')
%% Directory and filename
OriginalDirectory=cd;
WorkingDirectory = 'C:\'; 
cd(WorkingDirectory); 
  
%% How many dye kind??  Two channel 2 , three channel 3 %%
ColorNumber=2;  dye1c = 'g' ;    dye2c = 'r';  
  
%% Filenames     
filename_head = 'hel3dm'; 
    
%% Save Directgory  
SaveDirectory = 'C:\';
SaveDirectory = sprintf('%s\\%s', WorkingDirectory, filename_head);

%% Data Correcption %%    
dbackground_g=0; 
d2background_g=0; %-100;  
abackground_g=0; 
dbackground_r=0;  
d2background_r=0;  
abackground_r=0; 
leakage12=0.1;   %0.11 
leakage21=0; 
gamma12=1;  %1.0 

% leakage13=0.013; leakage 23=0.12; ga mma13=1; 

%% Trace Display Option 
BottomLimit_g=-100;   
UpperLimit_g=1000; 
BottomLimit_r=-100;
UpperLimit_r=1000;   
FirstNumber = 0; 
LastNumber = 10;     

BottomLimit_g=-500;  
UpperLimit_g=1500; 
BottomLimit_r=-500;
UpperLimit_r=1500;   
FirstNumber = 0; 
LastNumber = 10;     

BottomLimit_g=-500;  
UpperLimit_g=5500; 
BottomLimit_r=-500;
UpperLimit_r=5500;   
FirstNumber = 0; 
LastNumber = 10;     

%% Set Plotwidth %%   
PlotWidth = 0.15 %0.35 % standard = 0.8
PlotWidthStandard = 800;  
autoPlotWidth = 'y'
 
%% Draw Guide line?? %%
flow_frame1 = 30 / 1;  
flow_frame2 = flow_frame1; 
flow_frame3 = flow_frame1;

%% Movie required?? %%
DoseMovieNeed = 'y';
movie_initialframe = 3; 
%spotshape_x=[4.5 4.5 3.5 3.5 2.5 2.5 3.5 3.5 4.5 4.5 7.5 7.5 8.5 8.5 9.5 9.5 8.5 8.5 7.5 7.5 4.5]; % spotsize=7
%spotshape_y=[2.5 3.5 3.5 4.5 4.5 7.5 7.5 8.5 8.5 9.5 9.5 8.5 8.5 7.5 7.5 4.5 4.5 3.5 3.5 2.5 2.5]; % spotsize=7
spotshape_x=[4.5 4.5 3.5 3.5 4.5 4.5 7.5 7.5 8.5 8.5 7.5 7.5 4.5]; % spotsize=5
spotshape_y=[3.5 4.5 4.5 7.5 7.5 8.5 8.5 7.5 7.5 4.5 4.5 3.5 3.5];

%% Step Analysis %%
DoesStepAnalysisNeed = 'n';
ana_start = 25;
ana_end_set = 1000;
StepHeightMin = 100;

%% Filenames
filename_traces = [filename_head '.' num2str(ColorNumber) 'color_2alex_traces'];
filename_movie = [filename_head '.' num2str(ColorNumber) 'color_2alex_movies'];
filename_time = [filename_head '_time.dat'];
filename_add_nucleation = [ filename_head '_nucleation.dat']; 
filename_add = [ filename_head '_select.dat'];
filename_add_2regionA = [ filename_head '_2regionA.dat'];
filename_add_2regionB = [ filename_head '_2regionB.dat'];
filename_add_region = [ filename_head '_region.dat' ];
filename_add_region_data=[ filename_head '_region_data.dat' ];

%% binning required?? %%
DoseBinningNeed = 'n';
binwidth=5;

%% filter required?? %%
DoesFilterNeed = 'n';

%% Average and E level save?? %%
Is_Avg_and_E_save ='y';
r1 = 1;
r2 = 20;

%% Time unit is 'ms'? %%
Time_unit_ms = 'n';

%% Define time unit %%
fileinfo = dir([filename_head '.log']);
if sum(size(fileinfo)) == 1
    disp(['No log file : '  filename_head '.log']);
end

date = fileinfo.date;
fileid_log = fopen([filename_head '.log'],'r');		%% .log file
timeunit = textscan(fileid_log, '%*[Exposure time:] %f', 1);
timeunit = timeunit{1,1};
if Time_unit_ms =='y'
	timeunit= timeunit/1000;
end

% timeunit=0.05;     %% You want manual? Then use this line.
textscan(fileid_log, '%*[Acquisition mode:  Full 512x512 1x1 256x256 2x2 Binning]', 1);
gain = textscan(fileid_log, '%*[Gain:] %d', 1);
gain = gain{1,1};
scaler = textscan(fileid_log, '%*[Data scaler:] %f', 1);
scaler = scaler{1,1};
background = textscan(fileid_log, '%*[Background subtraction:] %d', 1);
background = background{1,1};
fclose(fileid_log);

%% Reading data %%
fileid = fopen(filename_traces,'r');	%% .trace file is made from idl(.run ana_all)
if fileid == -1
    disp(['No data file  '  filename_head]);
end

time_length = fread(fileid, 1, 'int32');
disp('The length of the time traces is: ')	% the length of the time trace
disp(time_length);
time = (0:(time_length-1))*timeunit;
if mod(time_length,2)==0
    time_g = (0:+2:(time_length-2))*timeunit;
    time_r = (1:+2:(time_length-1))*timeunit;
else
    time_g = (0:+2:(time_length-3))*timeunit;
    time_r = (1:+2:(time_length-2))*timeunit;
end

NumberofTraces = fread(fileid, 1, 'int16');
NumberofPeaks = NumberofTraces/ColorNumber;
disp('The number of traces and peaks are:')
disp(NumberofTraces);
disp(NumberofPeaks);
Data = fread(fileid, [NumberofTraces  time_length],'int16');
SpotDiameter = fread(fileid, 1, 'int16');
disp('Done reading trace data.');
fclose(fileid);

Temp_i = [];
Tempfirstpoint = [];
Templastpoint = [];
Templength = [];
Temptime_region = [];
TempFret12_region = [];
TempFret13_region = [];
TempFret23_region = [];
TempDonor_g_region = [];
TempDonor2_g_region = [];
TempAcceptor_g_region = [];
TempDonor_r_region = [];
TempDonor2_r_region = [];
TempAcceptor_r_region = [];
TempStepNumManual = [];
Temp_trace = [];
TempTime1 = [];
TempTime2 = [];
TempTime3 = [];
output = [];
output1 = [];
output2 = [];
output3 = [];
TempFret12a = [];
TempFret12b = [];        
TempDwellTime = [];     
TempAvgFRET = [];     
TempAvgSig1 = [];     
TempAvgSig2 = [];     
TempAvgSig3 = [];     
TempAvgSig4 = [];   
  
%% Movie Display
if DoseMovieNeed == 'y'
	cd(WorkingDirectory);
    disp(filename_movie);
	fileid_movie = fopen(filename_movie, 'r');
    if fileid_movie ~= -1
        peaks_total_width = fread(fileid_movie, 1, 'int16');
        peak_height = fread(fileid_movie, 1, 'int16');
        peaks_number = peaks_total_width/peak_height;
        file_information = dir(filename_movie);
        film_time_length = (file_information.bytes-4)/(peak_height*peaks_total_width);
        fclose(fileid_movie);

        disp('peaks_total_width, height, number, film_time_length: ');
        disp(peaks_total_width);
        disp(peak_height);
        disp(peaks_number);
        disp(film_time_length);
	
        %peak_line = zeros(ColorNumber*peak_height, 1, 'uint8');
        %%peaks=fread(fileid_movie, peaks_total_width * peak_height * film_time_length, 'uint8');
        peak = zeros(ColorNumber*peak_height, peak_height, 'uint8');

        if peaks_number ~= NumberofTraces
    		disp('error: Different trace numbers between .trace and .movies');
        	return;
        end

        if film_time_length ~= time_length
    		disp('error: Different time length between .trace and .movies');
        	return;
        end
    else
        DoseMovieNeed = 'n';
		disp('Movie viewer is turned off.');
    end
end


%% Convert raw data into donor and acceptor traces %%

if autoPlotWidth == 'y'
    PlotWidth = PlotWidth * (time_length / PlotWidthStandard);
    if PlotWidth>0.8
       PlotWidth=0.8 
    end
end

time_length_each = floor(time_length/2);
DonorRawData_1 = zeros(NumberofPeaks, time_length_each, 'double');
Donor2RawData_1 = zeros(NumberofPeaks, time_length_each, 'double');
AcceptorRawData_1 = zeros(NumberofPeaks, time_length_each, 'double');
DonorRawData_2 = zeros(NumberofPeaks, time_length_each, 'double');
Donor2RawData_2 = zeros(NumberofPeaks, time_length_each, 'double');
AcceptorRawData_2 = zeros(NumberofPeaks, time_length_each, 'double');
binlength = int32(time_length/binwidth-1);
bintime = zeros(binlength, 1, 'double');
binEraw = zeros(binlength, 1, 'double');
binEcorrect = zeros(binlength, 1, 'double');


for m=1:binlength
	bintime(m) = double(m-1)*(binwidth*timeunit);
end

% Rawdata를 나누어서 trace만들기
if ColorNumber == 2
    for i=1:NumberofPeaks
        for j=1:time_length_each
            DonorRawData_1(i,j) = Data(i*2-1,j*2-1);
            Donor2RawData_1(i,j) = Data(i*2,j*2-1);
            %AcceptorRawData_1(i,j) = Data(i*2,j*2-1);
            AcceptorRawData_1(i,j) = 0;
            DonorRawData_2(i,j) = Data(i*2-1,j*2);
            Donor2RawData_2(i,j) = Data(i*2,j*2);
            %AcceptorRawData_2(i,j) = Data(i*2,j*2);
            AcceptorRawData_2(i,j) = 0;
        end
    end    
end
if ColorNumber == 3
    for i=1:NumberofPeaks
        for j=1:time_length_each
            DonorRawData_1(i,j) = Data(i*3-2,j*2-1);
            Donor2RawData_1(i,j) = Data(i*3-1,j*2-1);
            AcceptorRawData_1(i,j) = Data(i*3,j*2-1);
            DonorRawData_2(i,j) = Data(i*3-2,j*2);
            Donor2RawData_2(i,j) = Data(i*3-1,j*2);
            AcceptorRawData_2(i,j) = Data(i*3,j*2);
        end
    end
end

clear Data;



%% calculate, plot and save average traces %%
%각 trace의 average를 구한후, 첫번째 laser excitation의 Cy3 channel과 두번째 laser
%excitation의 cy3 channel을 비교한 후, 큰 쪽의 trace를 green laser excitation으로 작은 쪽의
%trace를 red laser excitation으로 배정

%%%%%%%%%%%%% red-green excitate inversion
%%% Comparison %%%
if sum(sum(DonorRawData_1, 1)) > sum(sum(DonorRawData_2, 1))
    DonorRawData_g = DonorRawData_1;
    Donor2RawData_g = Donor2RawData_1;
    AcceptorRawData_g = AcceptorRawData_1;
    DonorRawData_r = DonorRawData_2; 
    Donor2RawData_r = Donor2RawData_2;
    AcceptorRawData_r = AcceptorRawData_2;
    color_order = 0;
else
    DonorRawData_g = DonorRawData_2;
    Donor2RawData_g = Donor2RawData_2;
    AcceptorRawData_g = AcceptorRawData_2;
    DonorRawData_r = DonorRawData_1;
    Donor2RawData_r = Donor2RawData_1;
    AcceptorRawData_r = AcceptorRawData_1;   
	color_order = 1;
end

DonorRawAverage_g = sum(DonorRawData_g, 1) / NumberofPeaks;
Donor2RawAverage_g = sum(Donor2RawData_g, 1) /NumberofPeaks;
AcceptorRawAverage_g = sum(AcceptorRawData_g, 1) / NumberofPeaks;
DonorRawAverage_r = sum(DonorRawData_r, 1) / NumberofPeaks;
Donor2RawAverage_r = sum(Donor2RawData_r, 1) /NumberofPeaks;
AcceptorRawAverage_r = sum(AcceptorRawData_r, 1) / NumberofPeaks;

clear DonorRawData_1;
clear Donor2RawData_1;
clear AcceptorRawData_1;
clear DonorRawData_2;
clear Donor2RawData_2;
clear AcceptorRawData_2;

scrsz = get(0,'ScreenSize');
figure('Name','Raw Data Ensenble Average','OuterPosition',[1 scrsz(4)/2 scrsz(4)/2 scrsz(4)/2]);
hdl1 = gcf;

%Green excitation의 signal average
subplot(2,1,1);
plot(time_g, DonorRawAverage_g - dbackground_g, dye1c, time_g, Donor2RawAverage_g - d2background_g, dye2c);
title(['Raw data Average signals in excitation laser ' dye1c]);
zoom on;

%Red excitation의 signal average
subplot(2,1,2);
plot(time_r, DonorRawAverage_r - dbackground_r, dye1c, time_r, Donor2RawAverage_r - d2background_r, dye2c);

title(['Raw data Average signals in excitation laser ' dye2c]);
zoom on;

if Is_Avg_and_E_save =='y'
    AverageOutput_g = [time_g' DonorRawAverage_g' Donor2RawAverage_g' AcceptorRawAverage_g'];
    AverageOutput_r = [time_r' DonorRawAverage_r' Donor2RawAverage_r' AcceptorRawAverage_r'];
    save([filename_head '_avg_g.dat'], 'AverageOutput_g', '-ascii');
    save([filename_head '_avg_r.dat'], 'AverageOutput_r', '-ascii');
    
    r1 = r1/2;
    r2 = r2/2;
    N_DonorRawAverage = DonorRawAverage_g / mean(DonorRawAverage_g(r1:r2));
    N_AcceptorRawAverage = Donor2RawAverage_g / mean(Donor2RawAverage_g(r1:r2));
    NormalizedAverageOutput = [time_g' N_DonorRawAverage' N_AcceptorRawAverage'];
    save([filename_head '_norm_avg_g.dat'], 'NormalizedAverageOutput', '-ascii'); 
    N_DonorRawAverage = DonorRawAverage_r / mean(DonorRawAverage_r(r1:r2));
    N_AcceptorRawAverage = Donor2RawAverage_r / mean(Donor2RawAverage_r(r1:r2));
    NormalizedAverageOutput = [time_r' N_DonorRawAverage' N_AcceptorRawAverage'];
    save([filename_head '_norm_avg_r.dat'], 'NormalizedAverageOutput', '-ascii');   
    
    % Intensity average of each molecule
    mol_DonorRawAverage_g = sum(DonorRawData_g, 2) / time_length_each;
    mol_Donor2RawAverage_g = sum(Donor2RawData_g, 2) /time_length_each;
    mol_AcceptorRawAverage_g = sum(AcceptorRawData_g, 2) / time_length_each;
    mol_DonorRawAverage_r = sum(DonorRawData_r, 2) / time_length_each;
    mol_Donor2RawAverage_r = sum(Donor2RawData_r, 2) /time_length_each;
    mol_AcceptorRawAverage_r = sum(AcceptorRawData_r, 2) / time_length_each;
    
    EachMoleculeAverageOutput = [mol_DonorRawAverage_g mol_Donor2RawAverage_g mol_DonorRawAverage_r mol_Donor2RawAverage_r];
    save(['avg_EachMolecule_' filename_head '.dat'], 'EachMoleculeAverageOutput', '-ascii');
    
end

%% calculate E level from the first 20 points and plot histograms of E level and total intensity. Also save the same info
tempDonor_g = reshape(DonorRawData_g(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
tempDonor2_g = reshape(Donor2RawData_g(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
tempAcceptor_g = reshape(AcceptorRawData_g(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
tempDonor_r = reshape(DonorRawData_r(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
tempDonor2_r = reshape(Donor2RawData_r(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
tempAcceptor_r = reshape(AcceptorRawData_r(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);

EachTotal_23_r = tempDonor2_r + tempAcceptor_r;
EachTotal_23_r = (EachTotal_23_r~=0).*EachTotal_23_r + (EachTotal_23_r==0)*1;	% remove zeros
EachTotal_123_g = tempDonor_g + tempDonor2_g + tempAcceptor_g;
EachTotal_123_g = (EachTotal_123_g~=0).*EachTotal_123_g + (EachTotal_123_g==0)*1;	% remove zeros

E_level_23 = tempAcceptor_r./EachTotal_23_r;
E_level_12 = tempDonor2_g./((1-E_level_23).*tempDonor_g + tempDonor2_g);
E_level_13 = (tempAcceptor_g - E_level_23.*(tempDonor2_g + tempAcceptor_g))./(tempDonor_g + tempAcceptor_g - E_level_23.*(EachTotal_123_g));

figure('Name','Raw data analysis: Early points','OuterPosition',[1 1 scrsz(4) scrsz(4)/2+100]);
hdl2 = gcf;

subplot(3,4,1); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_12,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
axis([-0.1 1.1 0 100]);
axis 'auto y'
title([ 'first ' num2str(FirstNumber) 'p Raw FRET12 histogram' ]);
zoom on;

subplot(3,4,2); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_13,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw FRET13 histogram' ]);
zoom on;

subplot(3,4,3); % 2*3 figure_ upper left (the last number shows the location of figure)
hist(E_level_23,-0.1:0.02:1.1); % histogram for first 10 point with the number of 50 shows the bin size
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw FRET23 histogram' ]);
zoom on;

subplot(3,4,4);
hist(EachTotal_123_g,-100:50:4000);
temp=axis;
temp(1)=-100;
temp(2)=4000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Total intensity histogram' ]);
zoom on;

subplot(3,4, [5 6 9 10]);
plot(E_level_12, EachTotal_123_g,'b+', 'MarkerSize', 2);
temp=axis;
temp(1)=-0.1;
temp(2)=1.1;
temp(3)=-100;
temp(4)=4000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Total Intensity vs. FRET' ]);
zoom on;

DonorFirstData_g = reshape(DonorRawData_g(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
Donor2FirstData_g = reshape(Donor2RawData_g(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
AcceptorFirstData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
size(DonorRawData_g)
time_length_each
DonorLastData_g = reshape(DonorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
Donor2LastData_g = reshape(Donor2RawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
AcceptorLastData_g = reshape(AcceptorRawData_g(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
DonorFirstData_r = reshape(DonorRawData_r(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
Donor2FirstData_r = reshape(Donor2RawData_r(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
AcceptorFirstData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, 2:FirstNumber+1), NumberofPeaks*FirstNumber, 1);
DonorLastData_r = reshape(DonorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
Donor2LastData_r = reshape(Donor2RawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);
AcceptorLastData_r = reshape(AcceptorRawData_r(1:NumberofPeaks, (time_length_each + 1 - LastNumber):time_length_each), NumberofPeaks*LastNumber, 1);


subplot(3,4,7);
hist(DonorFirstData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Donor Intensity histogram']);
zoom on;

subplot(3,4,8);
hist(DonorLastData_g,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Donor Intensity histogram']);
zoom on;

subplot(3,4,11);
hist(Donor2FirstData_r,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'first ' num2str(FirstNumber) 'p Raw Donor2 Intensity histogram']);
zoom on;

subplot(3,4,12);
hist(Donor2LastData_r,-300:5:2000);
temp=axis;
temp(1)=-300;
temp(2)=2000;
axis(temp);
title([ 'last ' num2str(LastNumber) 'p Raw Donor2 Intensity histogram']);
zoom on;

%% Start to servey
Temptime = [];
%TempFret = [];
%TempDonor = [];
%TempAcceptor = [];
TempFret12 = [];
TempFret13 = [];
TempFret23 = [];
TempDonor_g = [];
TempDonor2_g = [];
TempAcceptor_g = [];
TempDonor_r = [];
TempDonor2_r= [];
TempAcceptor_r = [];
output2 = [];
outputz = [];
outputz2 = [];

figure('Name','Trace analysis','OuterPosition',[scrsz(4)/2-200 0.05*scrsz(4) scrsz(4)+500 0.95*scrsz(4)]);
hdl_trace=gcf;

i=0;
history_n = 0;
history = zeros(1000, 1, 'int16');

while i < NumberofPeaks
    i = int32(int32(i) + 1);

    again=1;
	while ((again==1) && (i <= NumberofPeaks))
        if ColorNumber == 3
            DonorCorrect_g = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_g(i,:) - dbackground_g) - leakage21 * (Donor2RawData_g(i,:) - d2background_g));
            Donor2Correct_g = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_g(i,:) - d2background_g) - leakage12 * (DonorRawData_g(i,:) - dbackground_g)); 
    		AcceptorCorrect_g = gamma13 * ((AcceptorRawData_g(i,:) - abackground_g) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_g(i,:) - dbackground_g) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_g(i,:) -d2background_g));
			EachTotalCorrect_g = DonorCorrect_g + Donor2Correct_g + AcceptorCorrect_g;
			EachTotalCorrect_g = (EachTotalCorrect_g~=0).*EachTotalCorrect_g + (EachTotalCorrect_g==0)*0.0001;	% remove zeros
			
			DonorCorrect_r = ((1 + leakage12 + leakage13) / (1 - leakage12 * leakage21)) * ((DonorRawData_r(i,:) - dbackground_r) - leakage21 * (Donor2RawData_r(i,:) - d2background_r));
			Donor2Correct_r = gamma12 * ((1 + leakage21 + leakage23) / (1 - leakage12 * leakage21)) * ((Donor2RawData_r(i,:) - d2background_r) - leakage12 * (DonorRawData_r(i,:) - dbackground_r)); 
			AcceptorCorrect_r = gamma13 * ((AcceptorRawData_r(i,:) - abackground_r) - ((leakage23 * leakage12 - leakage13)/(1 - leakage12 * leakage21)) * (DonorRawData_r(i,:) - dbackground_r) + ((leakage13 * leakage21 - leakage23 ) / (1 - leakage12 * leakage21)) * (Donor2RawData_r(i,:) -d2background_r));
			AcceptorCorrect_r = AcceptorCorrect_r - direct * (Donor2Correct_r + AcceptorCorrect_r);
			EachTotalCorrect_r = Donor2Correct_r + AcceptorCorrect_r;
			EachTotalCorrect_r = (EachTotalCorrect_r~=0).*EachTotalCorrect_r + (EachTotalCorrect_r==0)*1;	% remove zeros

			Fret23 = AcceptorCorrect_r./EachTotalCorrect_r;
			Fret12 = Donor2Correct_g./((1-Fret23).*DonorCorrect_g + Donor2Correct_g);
			Fret13 = (AcceptorCorrect_g - Fret23.*(Donor2Correct_g + AcceptorCorrect_g))./(DonorCorrect_g + AcceptorCorrect_g - Fret23 .* (EachTotalCorrect_g));
        end
        if ColorNumber == 2
            DonorCorrect_g = (DonorRawData_g(i,:) - dbackground_g) + leakage12 * (DonorRawData_g(i,:) - dbackground_g);
            Donor2Correct_g = gamma12 * ((Donor2RawData_g(i,:) - d2background_g) - leakage12 * (DonorRawData_g(i,:) - dbackground_g));
            AcceptorCorrect_g = 0 * DonorRawData_g(i,:);
			EachTotalCorrect_g = DonorCorrect_g + Donor2Correct_g + AcceptorCorrect_g;
            EachTotalCorrect_g = (EachTotalCorrect_g~=0).*EachTotalCorrect_g + (EachTotalCorrect_g==0)*0.0001;	% remove zeros
			
            DonorCorrect_r = (DonorRawData_r(i,:) - dbackground_r) + leakage12 * (DonorRawData_r(i,:) - dbackground_r);
            Donor2Correct_r = gamma12 * ((Donor2RawData_r(i,:) - d2background_r) - leakage12 * (DonorRawData_r(i,:) - dbackground_r));
            AcceptorCorrect_r = 0 * DonorCorrect_r;
            EachTotalCorrect_r = Donor2Correct_r + AcceptorCorrect_r;
            EachTotalCorrect_r = (EachTotalCorrect_r~=0).*EachTotalCorrect_r + (EachTotalCorrect_r==0)*1;	% remove zeros
			
            Fret23 = 0 * DonorCorrect_g;
            Fret12 = Donor2Correct_g./EachTotalCorrect_g;
            Fret13 = 0 * DonorCorrect_g;
        end

        for j=2:time_length_each
            if (Fret13(j) < -0.3)
                Fret13(j) = 0;
            elseif (Fret13(j) > 1.1)
                Fret13(j) = 1;
            end
            if(Fret12(j) < -0.3)
                Fret12(j) = 0;
            elseif(Fret12(j) > 1.1)
                Fret12(j) = 1;
            end
            if(Fret23(j) < -0.3)
                Fret23(j) = 0;
            elseif(Fret23(j) > 1.1)
                Fret23(j) = 1;                
            end
        end
        
			again=0;

			EachTotalCorrect_g_normalized = EachTotalCorrect_g/max(EachTotalCorrect_g);
			DonorCorrect_g_Avg5point = smooth(DonorCorrect_g)';
			Donor2Correct_g_Avg5point = smooth(Donor2Correct_g)';
			Donor2Correct_r_Avg5point = smooth(Donor2Correct_r)';
		
	end

	%% Trace window
    figure(hdl_trace);
    
    %% trace1 Green laser excitation corrected trace
    subplot('position',[0.1 0.67 PlotWidth 0.30]); 
        
    plot(time_g, DonorCorrect_g, dye1c, time_g, Donor2Correct_g, dye2c); 
	if flow_frame1 ~= 0
        flow_time = flow_frame1*timeunit;
        line([flow_time flow_time],[BottomLimit_g UpperLimit_g],'Color','k', 'LineStyle', ':', 'LineWidth', 1);
    end
    
	if flow_frame2 ~= 0
        flow_time = flow_frame2*timeunit;
        line([flow_time flow_time],[BottomLimit_g UpperLimit_g],'Color','k', 'LineStyle', ':', 'LineWidth', 1);
    end
    
	if flow_frame3 ~= 0
        flow_time = flow_frame3*timeunit;
        line([flow_time flow_time],[BottomLimit_g UpperLimit_g],'Color','k', 'LineStyle', ':', 'LineWidth', 1);
    end
    
    temp=axis;
	temp(3)=BottomLimit_g;
	temp(4)=UpperLimit_g;
	grid on;
	axis(temp);
	title([dye1c ' Excitation      Molecule  ' num2str(i) '  / ' num2str(NumberofPeaks) '      file:' filename_head  '      Correction: ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(gamma12) '      ' date]);
    
	region_indice = find(Temp_i == i);
	number_region = size(region_indice);
	if number_region ~= 0
		FirstSelectX=[(Tempfirstpoint(region_indice)*timeunit)' (Tempfirstpoint(region_indice)*timeunit)'];
		FirstSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
		LastSelectX=[(Templastpoint(region_indice)*timeunit)' (Templastpoint(region_indice)*timeunit)'];
		LastSelectY=[(zeros(number_region) - 2)' (zeros(number_region) + 2)'];
    end

    %% trace2 Red laser excitation corrected trace
	
    subplot('position',[0.1 0.32 PlotWidth 0.30]); % 2nd
    plot(time_g, DonorCorrect_r, dye1c, time_g, Donor2Correct_r, dye2c); %, time_g, AcceptorCorrect_r, 'b');  
	if flow_frame1 ~= 0
        flow_time = flow_frame1*timeunit;
        line([flow_time flow_time],[BottomLimit_g UpperLimit_g],'Color','k', 'LineStyle', ':', 'LineWidth', 1);
    end
    
	if flow_frame2 ~= 0
        flow_time = flow_frame2*timeunit;
        line([flow_time flow_time],[BottomLimit_g UpperLimit_g],'Color','k', 'LineStyle', ':', 'LineWidth', 1);
    end
    
	if flow_frame3 ~= 0
        flow_time = flow_frame3*timeunit;
        line([flow_time flow_time],[BottomLimit_g UpperLimit_g],'Color','k', 'LineStyle', ':', 'LineWidth', 1);
    end
    
	temp=axis;
	temp(3)=BottomLimit_r;
	temp(4)=UpperLimit_r;
	grid on;
	axis(temp);
	title([dye2c ' Exciatation      Timeunit: ' num2str(timeunit) 'sec      Gain: ' num2str(gain) '      Scaler: ' num2str(scaler)  '      Spot diameter: ' num2str(SpotDiameter)]);

    region_indice = find(Temp_i == i);
	number_region = sum(size(region_indice)) - 1;
	FirstSelectX=[];
	FirstSelectY=[];
	LastSelectX=[];
	LastSelectY=[];
	for j=1:number_region
		tempx=[Tempfirstpoint(region_indice(j))*2*timeunit Tempfirstpoint(region_indice(j))*2*timeunit];
		tempy=[ (2 - 4 * mod(j,2)) (4 * mod(j,2) - 2)];
		FirstSelectX=[FirstSelectX tempx];
		FirstSelectY=[FirstSelectY tempy];
		tempx=[Templastpoint(region_indice(j))*2*timeunit Templastpoint(region_indice(j))*2*timeunit];
		LastSelectX=[LastSelectX tempx];
		LastSelectY=[LastSelectY tempy];
    end
    
    %% trace3 FRET trace
    
    subplot('position',[0.1 0.17 PlotWidth 0.10]);
    plot(time_g, Fret12, 'b');
	temp=axis;
	temp(3)=0;
	temp(4)=1;
	axis(temp);
    set(gca,'YTick', [0.2:0.2:1.2]);
    title('FRET');
	grid on;
	
    subplot('position',[0.93 0.42 0.05 0.15]);
	x = -0.1:0.02:1.1;
	[hX,hN]=hist(Fret12,x);
    barh(hN,hX,'k');
	temp=axis;
	temp(3)=-0.1;
	temp(4)=1.1;
	axis(temp);
	grid on;
	axis on;
   
   %% trace4
   
    subplot('position',[0.1 0.02 PlotWidth 0.10]);    

    GreenSmooth = smooth(DonorCorrect_g, 3);
    GreenTotal = DonorCorrect_g + Donor2Correct_g + 500;
    plot(time_g, DonorCorrect_g, dye1c, time_g, Donor2Correct_g, dye2c, time_g, GreenSmooth, 'k', time_g, GreenTotal, 'm'); % with smoothed
    temp=axis;
	temp(3)=BottomLimit_g;
	temp(4)=UpperLimit_g+500;
	grid on;
	axis(temp);
	title([dye1c ' Excitation      Molecule  ' num2str(i) '  / ' num2str(NumberofPeaks) '      file:' filename_head  '      Correction: ' num2str(dbackground_g) '  ' num2str(d2background_g) '  ' num2str(dbackground_r) '  ' num2str(d2background_r) '  ' num2str(leakage12) '  ' num2str(leakage21) '  ' num2str(gamma12) '      ' date]);

    %% Input Option
    
	if DoseMovieNeed == 'y'
        Xpoint = movie_initialframe; 
        
        fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
		for j=1:peak_height
			fseek(fileid_movie, startpoint , 'bof');
			peak_line=fread(fileid_movie, peak_height*3, 'uint8');
			peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
			startpoint = startpoint + peaks_total_width;
		end
		fclose(fileid_movie);
		
		subplot('position',[0.93 0.80 0.05 0.18]);
		colormap(hot);
		image(peak);
		axis off;
        hold on;
        plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
        hold off; 
        title([num2str(Xpoint) 'fr (' num2str(Xpoint*timeunit) ' sec)']);
        
		fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
        startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
		for j=1:peak_height
			fseek(fileid_movie, startpoint , 'bof');
            peak_line=fread(fileid_movie, peak_height*3, 'uint8');
			peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
			startpoint = startpoint + peaks_total_width;
		end
		fclose(fileid_movie);

		subplot('position',[0.93 0.60 0.05 0.18]);
		colormap(hot);
		image(peak);
		axis off;
        hold on;
        plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
        hold off; 
    end
    
    again=1;
	while again==1
		again=0;
		disp([num2str(i) '(' num2str( ceil( cast(i,'double') / NumberofPeaks * 10^4) / 10^2 ) '%)' ' (s=save by png, p=display movie)']);
		keyanswer =input('(t=terminate program, b=back, g=go) : ','s');
		answer = sscanf(keyanswer, '%s %*s');
		numberofanswer = sscanf(keyanswer, '%*s %f');
        
        if answer=='z'  %save a trace with a log for further analysis
            filename_add_selected = ['sel_AlexTrace_', filename_head, '_', num2str(i), '.dat'];
            outputz = [time_g' Fret12' DonorCorrect_g' Donor2Correct_g' DonorCorrect_r' Donor2Correct_r'];
            save(filename_add_selected, 'outputz', '-ascii');            
        end
        
        if answer=='@'  %reassign background
            BottomLimit_g=str2double(input('BottomLimit_g? ','s'));
            UpperLimit_g=str2double(input('UpperLimit_g? ','s'));
            BottomLimit_r=str2double(input('BottomLimit_r? ','s'));
            UpperLimit_r=str2double(input('UpperLimit_r? ','s'));
        end
        
		if answer=='o'
			again=1;
			[Xc,Yc] = ginput(1);
            Xc(1); Yc(1);
			Xpoint = round(Xc(1)/timeunit);

           	if DoseMovieNeed == 'y'
				fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
				
				startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
				for j=1:peak_height
					fseek(fileid_movie, startpoint, 'bof');
					peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
					peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
					startpoint = startpoint + peaks_total_width;
				end
				fclose(fileid_movie);
				
                subplot('position',[0.93 0.80 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                hold on;
                plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                hold off;  
                title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
        
				fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*3, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
				end
				fclose(fileid_movie);
				
                subplot('position',[0.93 0.60 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                hold on;
                plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                hold off;  
            end
        end
        
		if answer=='i'
%			again=1;
			[Xc,Yc] = ginput(1);
            Xc(1); Yc(1);
			Xpoint = round(Xc(1)/timeunit);
            
            if Xpoint>=0
                if DoseMovieNeed == 'y'
                    fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');

                    startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
                    for j=1:peak_height
                        fseek(fileid_movie, startpoint, 'bof');
                        peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                        peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                        startpoint = startpoint + peaks_total_width;
                    end
                    fclose(fileid_movie);

                    subplot('position',[0.93 0.80 0.05 0.18]);
                    colormap(hot);
                    image(peak);
                    axis off;
                    hold on;
                    plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                    hold off;  
                    title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);

                    fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
                    startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
                    for j=1:peak_height
                        fseek(fileid_movie, startpoint , 'bof');
                        peak_line=fread(fileid_movie, peak_height*3, 'uint8');
                        peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                        startpoint = startpoint + peaks_total_width;
                    end
                    fclose(fileid_movie);

                    subplot('position',[0.93 0.60 0.05 0.18]);
                    colormap(hot);
                    image(peak);
                    axis off;
                    hold on;
                    plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                    hold off;  
                end

                filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_g_', num2str(Xpoint,'%03.f'), '.png'];
                saveas(gcf, filename_add_selected)
            end
        end        
        
if answer=='p' 
			again=1;
			[Xc,Yc] = ginput(1);
            Xc(1); Yc(1);
			Xpoint = round(Xc(1)/timeunit);
            
            if Xpoint>=0
                if DoseMovieNeed == 'y'
                    fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');

                    startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
                    for j=1:peak_height
                        fseek(fileid_movie, startpoint, 'bof');
                        peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
                        peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                        startpoint = startpoint + peaks_total_width;
                    end
                    fclose(fileid_movie);

                    subplot('position',[0.93 0.80 0.05 0.18]);
                    colormap(hot);
                    image(peak);
                    axis off;
                    hold on;
                    plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                    hold off;  
                    title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);

                    fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
                    startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
                    for j=1:peak_height
                        fseek(fileid_movie, startpoint , 'bof');
                        peak_line=fread(fileid_movie, peak_height*3, 'uint8');
                        peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                        startpoint = startpoint + peaks_total_width;
                    end
                    fclose(fileid_movie);

                    subplot('position',[0.93 0.60 0.05 0.18]);
                    colormap(hot);
                    image(peak);
                    axis off;
                    hold on;
                    plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                    hold off;  
                end

                keyanswer=input('Save? ','s');
                answer_for_save = sscanf(keyanswer, '%s %*s');
                
                if Xpoint < 30
             		disp('under 30 Error !!!');
             		disp('under 30 Error !!!');
             		disp('under 30 Error !!!');
               end

                if answer_for_save == 'p'
                elseif answer_for_save == '1'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_in_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                elseif answer_for_save == 'q'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_pife_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                    elseif answer_for_save == 'w'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_tm_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                    elseif answer_for_save == 'e' 
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_rt_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                else
            		disp('Error !!!');
            		disp('Error !!!');
            		disp('Error !!!');
%                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_', answer_for_save, '_', num2str(Xpoint,'%03.f'), '.png'];
%                    saveas(gcf, filename_add_selected)
                end
            end
        end
        
        if answer=='a' 
                filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_PIFE.png'];
                saveas(gcf, filename_add_selected)
        end        
        
        if answer=='r'
			again=1;

            if DoseMovieNeed == 'y'
				fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
				
				startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + color_order*peak_height*peaks_total_width;
				for j=1:peak_height
					fseek(fileid_movie, startpoint, 'bof');
					peak_line=fread(fileid_movie, peak_height*ColorNumber, 'uint8');
					peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
					startpoint = startpoint + peaks_total_width;
				end
				fclose(fileid_movie);
				
                subplot('position',[0.93 0.80 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                hold on;
                plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                hold off;  
                title(['frame: ' num2str(Xpoint) '  time ' num2str(Xpoint*timeunit)]);
        
				fileid_movie = fopen([filename_head '.' num2str(ColorNumber) 'color_2alex_movies'], 'r');
                startpoint = 4 + uint32((i-1)*ColorNumber*peak_height) + 2*(uint32(Xpoint/2)-1)*peak_height*peaks_total_width + (1-color_order)*peak_height*peaks_total_width;
                for j=1:peak_height
                    fseek(fileid_movie, startpoint , 'bof');
                    peak_line=fread(fileid_movie, peak_height*3, 'uint8');
                    peak(1:peak_height*ColorNumber, j) = peak_line(1:peak_height*ColorNumber);
                    startpoint = startpoint + peaks_total_width;
				end
				fclose(fileid_movie);
				
                subplot('position',[0.93 0.60 0.05 0.18]);
                colormap(hot);
                image(peak);
                axis off;
                hold on;
                plot(spotshape_x, spotshape_y, spotshape_x, spotshape_y+11, 'color', 'g');
                hold off;  
            end
            
            keyanswer=input('Save? ','s');
        	answer_for_save = sscanf(keyanswer, '%s %*s');
        
                if answer_for_save == 'p'
                elseif answer_for_save == '1'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_in_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                elseif answer_for_save == 'q'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_pife_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                    elseif answer_for_save == 'w'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_tm_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                    elseif answer_for_save == 'e'
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_rt_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                else
                    filename_add_selected = [SaveDirectory, '\', num2str(i,'%03.f'), '_', answer_for_save, '_', num2str(Xpoint,'%03.f'), '.png'];
                    saveas(gcf, filename_add_selected)
                end
        end        
        disp('again end');
    end
        
	if answer=='b'
		if i>1 && history_n > 0
			while history(history_n)==i
				history_n = history_n - 1;
			end
			i=history(history_n)-1;
			history_n = history_n - 1;
		else
			i=i-1;
		end
	else
		history_n = history_n + 1;
		history(history_n)=i;
    end
    
	if answer=='g'
		answer=input('number to go : ','s');
		gonumber = str2num(answer);
		if gonumber > 0 && gonumber <= NumberofPeaks
            i = gonumber - 1;
		end
	end

	if answer=='t'
		close all;
		cd(OriginalDirectory);
		break;
	end
end

disp('Program end.')

cd(OriginalDirectory);

clear all;
close all;
end
