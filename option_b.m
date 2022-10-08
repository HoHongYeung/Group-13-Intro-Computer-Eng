clear;
fid = fopen('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\listData.txt','rt');
textLine = fscanf(fid,'%s',1);
ind_s = 0;
ind_aa = 0;
num_s=0;
num_aa=0;
Orig_s=[];
Orig_aa=[];

while (strcmp(textLine,'.')~=1)
    filename=strcat('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\labels\',textLine,'.lab');
    file=fopen(filename,'rt');
    textCol = textscan(file,'%f %f %s');
    audioname=strcat('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavOrig\',textLine,'.wav');
    [inpSigWav,Fs] = audioread(audioname);
    
    
    % Store the data into a cell/standard array
    for i=1:size(textCol{1},1)
        if strcmp(textCol{3}(i),'s')==1
            ind_s=ind_s+1;
            num_s=num_s+1;
            timePhStart_s(ind_s) = textCol{1}(i);
            timePhEnd_s(ind_s) = textCol{2}(i);
            Orig_s{num_s}=inpSigWav;
        end
        
        if strcmp(textCol{3}(i),'aa')==1
            ind_aa=ind_aa+1;
            num_aa=num_aa+1;
            timePhStart_aa(ind_aa) = textCol{1}(i);
            timePhEnd_aa(ind_aa) = textCol{2}(i);
            Orig_aa{num_aa}=inpSigWav;
        end
    end
    
        
        
    fclose(file);
    textLine = fscanf(fid,'%s',1);
end

timeSegStart_ms_s = timePhStart_s*10^(-4) + (timePhEnd_s*10^(-4) - timePhStart_s*10^(-4))/2 - 12.5;
timeSegEnd_ms_s = timeSegStart_ms_s + 25;
indexSegStart_s = round(timeSegStart_ms_s*Fs/1000,0);
indexSegEnd_s = round(timeSegEnd_ms_s*Fs/1000,0);

timeSegStart_ms_aa = timePhStart_aa*10^(-4) + (timePhEnd_aa*10^(-4) - timePhStart_aa*10^(-4))/2 - 12.5;
timeSegEnd_ms_aa = timeSegStart_ms_aa + 25;
indexSegStart_aa = round(timeSegStart_ms_aa*Fs/1000,0);
indexSegEnd_aa = round(timeSegEnd_ms_aa*Fs/1000,0);

for j=1:ind_s
    segOrig_phS(j,1:201)=(Orig_s{j}(indexSegStart_s(j):indexSegEnd_s(j),1))';
end

for k=1:ind_aa
    segOrig_phAA(k,1:201)=(Orig_aa{k}(indexSegStart_aa(k):indexSegEnd_aa(k),1))';
end

save('segAllData.mat','segOrig_phAA','segOrig_phS');
load ('segAllData.mat');

fclose(fid);

%Plot the figure

[inpSigWav1,Fs] = audioread('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavOrig\MDPK0\SA1.wav');
fid1 = fopen('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\listData.txt','rt');
textLine1 = fscanf(fid,'%s',1);
ind_s_sa1=0;
ind_aa_sa1=0;

while (strcmp(textLine1,'.')~=1)
    filename1=strcat('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\labels\',textLine1,'.lab');
    file1=fopen(filename1,'rt');
    textCol1 = textscan(file1,'%f %f %s');
    
    % Store the data into a cell/standard array
    for m=1:size(textCol1{1},1)
        if strcmp(textCol1{3}(m),'s')==1
            ind_s_sa1=ind_s_sa1+1;
            if (strcmp(textLine1,'MDPK0\SA1')==1)
                signal_sa1_s=segOrig_phS(ind_s_sa1,1:end);
                x1_s=indexSegStart_s(ind_s_sa1);
                x2_s=indexSegEnd_s(ind_s_sa1);
                break;
            end
        end
    end
    
    for n=1:size(textCol1{1},1)
        if strcmp(textCol1{3}(n),'aa')==1
            ind_aa_sa1=ind_aa_sa1+1;
            if (strcmp(textLine1,'MDPK0\SA1')==1)
                signal_sa1_aa=segOrig_phAA(ind_aa_sa1,1:end);
                x1_aa=indexSegStart_aa(ind_aa_sa1);
                x2_aa=indexSegEnd_aa(ind_aa_sa1);
                break;
            end
        end
    end
    fclose(file1);
    textLine1 = fscanf(fid1,'%s',1);
end

fclose(fid1);

figure;
plot([8189:8389],signal_sa1_s(1:201),'r');
grid; xlabel('Frequncy Index'); ylabel('Amplitude');
legend('signal segment of s');


figure;
plot([6436:6636],signal_sa1_aa(1:201),'g');
grid; xlabel('Frequncy Index'); ylabel('Amplitude');
legend('signal segment of aa');








