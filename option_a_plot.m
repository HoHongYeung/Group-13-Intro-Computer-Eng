clear;
fileName = 'C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavOrig\MCPM0\sa1.wav';
[inpSigWav,Fs] = audioread(fileName);%  Load the input signal
nSamples = length(inpSigWav);
    
%  Define the FIR filter coefficients
firCoef = [-0.08,0.24,0.4,0.4,0.16,-0.24,0.08];
impres=[-0.08,0.24,0.4,0.4,0.16,-0.24,0.08,zeros(1, 249)];
nCoef = length(firCoef);
IMPRES=abs(fft(impres));
magIMPRES=20 * log10(IMPRES);
    
%  Add zero samples to the beginning and to the end of the input signal to take account of the filter length
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
plot([0:127]*8000/256,magIMPRES(1:128));
grid; xlabel('Frequncy'); ylabel('Amplitude[db]');
figure;
plot([1900 : 2300], real(inpSigWav(1900:2300)),'r',[1900 : 2300],real(outSigWav(1900:2300)),'g');
% plot([1900:2300],real(20*log10(INPSIGWAV(1900:2300))),'g',[1900:2300],real(20*log10(OUTSIGWAV(1900:2300))),'b'); 
grid; xlabel('Frequncy Index'); ylabel('Amplitude');
legend('Input signal', 'Output signal')
    
%  Store the output signal into .wav file
fileNameOut =  'C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\outputSignalatest2MCPM0sa1.wav';
audiowrite(fileNameOut, outSigWav, Fs);