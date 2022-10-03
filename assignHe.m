clear;
fid = fopen('listData.txt','rt');
ind=0;

% Read the column(column 1)
textCol1 = fscanf(fid,'%s',1);
while (strcmp(textCol1,'.')~=1)   %While the content of column 1 is not '.'
    textLine1 = fscanf(fid,'%s',1);

    ind = ind+1;

    content1=textLine1(1:6);     %分段处理，头部不变
    content2=textLine1(7:end);   %分段处理，音频文件名
    
    Enstr=strcat(content2,'.wav');
    fileName=strcat('wavOrig','\',content1,'\',Enstr);

    textLine1(7:8)=lower(textLine1(7:8)); %音频文件名小写



  [inpSigWav,Fs] = audioread(fileName);  %读取文件
  nSamples = length(inpSigWav);

  firCoef = [-0.08, 0.24,0.4,0.4,0.16,-0.24,0.08];
  nCoef = length(firCoef);

    inpSigWavExt = [zeros(1,1); inpSigWav; zeros(1,1)];
    nSamplesNew = length(inpSigWavExt);
    %  Allocate memory for the output signal
    outSigWav = zeros(nSamplesNew,1);
    %  Perform filtering
    for i=nCoef:nSamplesNew,
        outSigWav(i) = firCoef(1)*inpSigWavExt(i)+firCoef(2)*inpSigWavExt(i-1)+firCoef(3)*inpSigWavExt(i-2)+firCoef(4)*inpSigWavExt(i-3)+firCoef(5)*inpSigWavExt(i-4)+firCoef(6)*inpSigWavExt(i-5)+firCoef(7)*inpSigWavExt(i-6);
    end
    INPSIGWAV=fft(inpSigWav);
    OUTSIGWAV=fft(outSigWav);
    %  Plot the input and output signals
    figure;
    plot([1900:2300],real(20*log10(INPSIGWAV(1900:2300))),'g',[1900:2300],real(20*log10(OUTSIGWAV(1900:2300))),'b'); 
    grid; xlabel('Index'); ylabel('Amplitude[db]');
    legend('Input signal', 'Output signal')
    
    %  Store the output signal into .wav file
    fileNameOut = strcat('wavOrig','\',content1,'\',Enstr);
    audiowrite(fileNameOut, outSigWav, Fs);
end
fclose(fid);


A题ALL
clear;
fid = fopen('listData.txt','rt');
ind=0;

% Read the column(column 1)
textCol1 = fscanf(fid,'%s',1);
while (strcmp(textCol1,'.')~=1)  
    textLine1 = fscanf(fid,'%s',1);

    ind = ind+1;
    cont='textLine1(1:6).textLine1(7:8).textLine1(9:end)';
     endtr=textCol1(7:8);
    if endtr>='A' & endtr<='Z'
        lcl=lower(endtr);
    end
    cell_str=strsplit(cont,'.');
    path='E:\matlab\LMexercise\wavOrig\';
    fileName=strcat(save_path,textCol1(1:6),endtr,textCol1(9:end),'.wav');

    storeDataStruct(ind).name = fileName;

    inpSigWav,Fs] = audioread(fileName);  %读取文件
  nSamples = length(inpSigWav);

   firCoef = [-0.08, 0.24,0.4,0.4,0.16,-0.24,0.08];
   longfirCoef=[zeros(1,125),firCoef, zeros(1,125)];
   nCoef = length(firCoef);
  

    inpSigWavExt = [zeros(1,1); inpSigWav; zeros(1,1)];
    nSamplesNew = length(inpSigWavExt);
    %  Allocate memory for the output signal
    outSigWav = zeros(nSamplesNew,1);
    %  Perform filtering
    for i=nCoef:nSamplesNew,
        outSigWav(i) = firCoef(1)*inpSigWavExt(i)+firCoef(2)*inpSigWavExt(i-1)+firCoef(3)*inpSigWavExt(i-2)+firCoef(4)*inpSigWavExt(i-3)+firCoef(5)*inpSigWavExt(i-4)+firCoef(6)*inpSigWavExt(i-5)+firCoef(7)*inpSigWavExt(i-6);
    end

    %  Plot the input and output signals
    figure;
    plot([1900:2300],real(20*log10(inpSigWav(1900:2300))),'g',[1900:2300],real(20*log10(outSigWav(1900:2300))),'b'); 
    grid; xlabel('Index'); ylabel('Amplitude[db]');
    legend('Input signal', 'Output signal')
    
    %  Store the output signal into .wav file
      fileNameOut = strcat('wavOrig','\',cont,'\',Enstr,'.wav');
      audiowrite(fileNameOut, outSigWav, Fs);
end
    % fir
    figure;
    magFreqChar=abs(fft(longfirCoef));
    magFreqChar_dB=20*log10(magFreqChar);
    plot(magFreqChar_dB);
    plot([0:128]*8000/257,magFreqChar_dB(1:129));

fclose(fid);

Another way：
clear;
fid = fopen('listData.txt','rt');
ind=0;

% disp(textLine1);
% syms textLine1;

% Read the column(column 1)
textCol1 = fscanf(fid,'%s',1);
while (strcmp(textCol1,'.')~=1)   %While the content of column 1 is not '.'
    
    ind = ind+1;
    save_path='E:\matlab\LMexercise\wavOrig\';
    content='textLine1(1:6).textLine1(7:8).textLine1(9:end)';
    cell_str=strsplit(content,'.');
    endtr=textCol1(7:8);
    if endtr>='A' & endtr<='Z'
        endtrs=lower(endtr);
    end
        
    sep_name=cell_str{1,1};
    newname=strcat(save_path,textCol1(1:6),endtrs,textCol1(9:end),'.wav');
    % Store the data into a structure array
    storeDataStruct(ind).name = newname;
    %     textCol1 = fscanf(fid,'%s',1);
    
    [inpSigWav,Fs] = audioread(newname);%  Load the input signal
    nSamples = length(inpSigWav);
    
    %  Define the FIR filter coefficients
    firCoef = [-0.08,0.24,0.4,0.4,0.16,-0.24,0.08];
    nCoef = length(firCoef);
    
    %  Add zero samples to the beginning and to the end of the input signal to take account of the filter length
    inpSigWavExt = [zeros(1,1); inpSigWav; zeros(1,1)];
    nSamplesNew = length(inpSigWavExt);
    %  Allocate memory for the output signal
    outSigWav = zeros(nSamplesNew,1);
    %  Perform filtering
    for i=nCoef:nSamplesNew,
        outSigWav(i) = firCoef(1)*inpSigWavExt(i)+firCoef(2)*inpSigWavExt(i-1)+firCoef(3)*inpSigWavExt(i-2)+firCoef(4)*inpSigWavExt(i-3)+firCoef(5)*inpSigWavExt(i-4)+firCoef(6)*inpSigWavExt(i-5)+firCoef(7)*inpSigWavExt(i-6);
    end

    
    %  Store the output signal into .wav file
    fileNameOut = strcat('E:\matlab\LMexercise\wavFilt\',textCol1(1:6),endtr,textCol1(9:end),'.wav');
    audiowrite(fileNameOut, outSigWav, Fs);
    
    % Read column 1 - next row
    textCol1 = fscanf(fid,'%s',1);
    
end

%  Define the FIR filter coefficients
firCoef = [-0.08,0.24,0.4,0.4,0.16,-0.24,0.08];
impres=[-0.08,0.24,0.4,0.4,0.16,-0.24,0.08,zeros(1, 249)];
nCoef = length(firCoef);
IMPRES=abs(fft(impres));
magIMPRES=20 * log10(IMPRES);
figure;
plot([0:127]*8000/256,magIMPRES(1:128));
grid; xlabel('Frequncy'); ylabel('Amplitude[db]');

fileNameOrig = 'E:\matlab\LMexercise\wavOrig\MCPM0\sa1.wav';
[inpSigWav,Fs] = audioread(fileNameOrig);%  Load the input signal

fileNameFilt = 'E:\matlab\LMexercise\wavFilt\MCPM0\sa1.wav';
[outSigWav,Fs] = audioread(fileNameFilt);%  Load the input signal

figure;
plot([1900 : 2300], real(inpSigWav(1900:2300)),'r',[1900 : 2300],real(outSigWav(1900:2300)),'g');
grid; xlabel('Index'); ylabel('Amplitude');
legend('Input signal', 'Output signal')

fclose(fid);
