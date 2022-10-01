clc
clear all
textLine1='r';
fid = fopen('C:\Users\SHLT\Desktop\listData.txt','rt');
while (strcmp(textLine1,'.')~=1)
  textLine1 = fscanf(fid,'%s',1);
  textLine1(7)=lower(textLine1(7));
  textLine1(8)=lower(textLine1(8));
  aa=textLine1(1:5);
  bb=textLine1(7:end);
  bb=strcat(bb,'.wav');
  fileName=strcat('C:\Users\SHLT\Desktop\wavOrig','\',aa,'\',bb);
  [inpSigWav,Fs] = audioread(fileName);
  nSamples = length(inpSigWav);

  firCoef = [-0.08, 0.24,0.4,0.4,0.16,-0.24,0.08];
  nCoef = length(firCoef);

  %  Add zero samples to the beginning and to the end of the input signal to 
  %  take account of the filter length
  inpSigWavExt = [zeros(nCoef-1,1); inpSigWav; zeros(nCoef-1,1)];
  nSamplesNew = length(inpSigWavExt);
  outSigWav = zeros(nSamplesNew,1);
  for i=nCoef:nSamplesNew
    outSigWav(i) = firCoef(1)*inpSigWavExt(i)+firCoef(2)*inpSigWavExt(i-1)+firCoef(3)*inpSigWavExt(i-2)+firCoef(4)*inpSigWavExt(i-3)+firCoef(5)*inpSigWavExt(i-4)+firCoef(6)*inpSigWavExt(i-5);
  end
  %定义滤波器的脉冲响应
  figure;
  plot([1001:1200],inpSigWav(1001:1200),'g',[1001:1200],outSigWav(1001:1200),'b'); 
  grid; xlabel('Sample index'); ylabel('Amplitude');
  legend('Input signal', 'Output signal')
%   bb(1)=upper(bb(1));
%   bb(2)=upper(bb(2));
  fileNameOut = strcat('C:\Users\SHLT\Desktop\wavFilt','\',aa,'\',bb);
  audiowrite(fileNameOut, outSigWav, Fs);
end