clear;
fid = fopen('listData.txt','rt');
ind=0;

% disp(textLine1);
% syms textLine1;

% Read the column(column 1)
textCol1 = fscanf(fid,'%s',1);
while (strcmp(textCol1,'.')~=1)%While the content of column 1 is not '.'
    textLine1 = fscanf(fid,'%s',1);

    ind = ind+1;
    save_path='C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavOrig\';
    content='textLine1(1:6).textLine1(7:8).textLine1(9:end)';
    cell_str=strsplit(content,'.');
    ucl=textLine1(7:8);
    if ucl>='A' & ucl<='Z'
        lcl=lower(ucl);
    end
        
    sep_name=cell_str{1,1};
    newname=strcat(save_path,textLine1(1:6),lcl,textLine1(9:end),'.wav');
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
    INPSIGWAV=fft(inpSigWav);
    OUTSIGWAV=fft(outSigWav);
    %  Plot the input and output signals
    figure;
    plot([700:900],real(20*log10(INPSIGWAV(700:900))),'g',[700:900],real(20*log10(OUTSIGWAV(700:900))),'b'); 
    grid; xlabel('Index'); ylabel('Amplitude[db]');
    legend('Input signal', 'Output signal')
    
    %  Store the output signal into .wav file
    fileNameOut = strcat('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavFilt\',textLine1(1:6),lcl,textLine1(9:end),'.wav');
    audiowrite(fileNameOut, outSigWav, Fs);
end
fclose(fid);