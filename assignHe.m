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
