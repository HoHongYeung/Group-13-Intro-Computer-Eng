clear;
fid = fopen('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\listData.txt','rt');
textLine = fscanf(fid,'%s',1);
ind_s = 0;
ind_aa = 0;
num_s=0;
num_aa=0;
Orig_s=[];
Orig_aa=[];%Store full audio

while (strcmp(textLine,'.')~=1)
    filename=strcat('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\labels\',textLine,'.lab');
    file=fopen(filename,'rt');
    textCol = textscan(file,'%f %f %s');%Scan the 3 columns
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

ind1_A=ceil(201*400/Fs);
ind2_A=ceil(201*1600/Fs);
ind1_B=ceil(201*2400/Fs);
ind2_B=ceil(201*4000/Fs);

for j=1:ind_s
    value1(j,1)=0;
    value1(j,2)=0;%for frequency region B
    %intercept the Orig_s according to the time index, and store 
    %the result in segOrig_phS
    segOrig_phS(j,1:201)=(Orig_s{j}(indexSegStart_s(j):indexSegEnd_s(j),1))';
    segOrig_phS_ft(j,1:201)=abs(fft(segOrig_phS(j,1:201)));
    %the |X(k)|^2 of the number of columns 11 to 41 of 
    %each row is added and stored into the first column of value1
    for m=ind1_A:ind2_A
        value1(j,1)=abs(segOrig_phS_ft(j,m))^2+value1(j,1);
    end
    aveEn_dB_A_s(j,1)=10*log10(value1(j,1)/(ind2_A-ind1_A+1));
    enOrigRegAB_phS(j,1)=aveEn_dB_A_s(j,1);
    for n=ind1_B:ind2_B
        value1(j,2)=abs(segOrig_phS_ft(j,n))^2+value1(j,2);
    end
    aveEn_dB_B_s(j,1)=10*log10(value1(j,2)/(ind2_B-ind1_B+1));
    enOrigRegAB_phS(j,2)=aveEn_dB_B_s(j,1);
end

for k=1:ind_aa
    value2(k,1)=0;
    value2(k,2)=0;
    segOrig_phAA(k,1:201)=(Orig_aa{k}(indexSegStart_aa(k):indexSegEnd_aa(k),1))';
    segOrig_phAA(k,1:201)=(Orig_aa{k}(indexSegStart_aa(k):indexSegEnd_aa(k),1))';
    segOrig_phAA_ft(k,1:201)=abs(fft(segOrig_phAA(k,1:201)));
    for o=ind1_A:ind2_A
        value2(k,1)=abs(segOrig_phAA_ft(k,o))^2+value2(k,1);
    end
    aveEn_dB_A_aa(k,1)=10*log10(value2(k,1)/(ind2_A-ind1_A+1));
    enOrigRegAB_phAA(k,1)=aveEn_dB_A_aa(k,1);
    for p=ind1_B:ind2_B
        value2(k,2)=abs(segOrig_phAA_ft(k,p))^2+value2(k,2);
    end
    aveEn_dB_B_aa(k,1)=10*log10(value2(k,2)/(ind2_B-ind1_B+1));
    enOrigRegAB_phAA(k,2)=aveEn_dB_B_aa(k,1);
end

fclose(fid);

sum_element_A_s=0;
diff_element_A_s=0;
for q=1:80% Calculate the mean value
    sum_element_A_s=aveEn_dB_A_s(q)+sum_element_A_s;
    m_A_s=sum_element_A_s/80;
end
for q=1:80% Calculating the variance and standard deviation
    diff_element_A_s=(aveEn_dB_A_s(q)-m_A_s)^2+diff_element_A_s;
    v_A_s=diff_element_A_s/80;
    sv_A_s=v_A_s^(0.5);
end
for q=1:80% Calculated by the formula of Gaussian PDF, and the results are 
          % stored in the first column of gPDF_s in the 159th line, 
          % merging all the results into a .mat file in the 170th line.
    gPDF_A_s(q)=(1/sqrt(2*pi*v_A_s))*exp(-(aveEn_dB_A_s(q)-m_A_s)^2/2/v_A_s);
end

sum_element_B_s=0;
diff_element_B_s=0;
for q=1:80
    sum_element_B_s=aveEn_dB_B_s(q)+sum_element_B_s;
    m_B_s=sum_element_B_s/80;
end
for q=1:80
    diff_element_B_s=(aveEn_dB_B_s(q)-m_B_s)^2+diff_element_B_s;
    v_B_s=diff_element_B_s/80;
    sv_B_s=v_B_s^(0.5);
end
for q=1:80
    gPDF_B_s(q)=(1/sqrt(2*pi*v_B_s))*exp(-(aveEn_dB_B_s(q)-m_B_s)^2/2/v_B_s);
end

sum_element_A_aa=0;
diff_element_A_aa=0;
for r=1:64
    sum_element_A_aa=aveEn_dB_A_aa(r)+sum_element_A_aa;
    m_A_aa=sum_element_A_aa/64;
end
for r=1:64
    diff_element_A_aa=(aveEn_dB_A_aa(r)-m_A_aa)^2+diff_element_A_aa;
    v_A_aa=diff_element_A_aa/64;
    sv_A_aa=v_A_aa^(0.5);
end
for r=1:64
    gPDF_A_aa(r)=(1/sqrt(2*pi*v_A_aa))*exp(-(aveEn_dB_A_aa(r)-m_A_aa)^2/2/v_A_aa);
end

sum_element_B_aa=0;
diff_element_B_aa=0;
for r=1:64
    sum_element_B_aa=aveEn_dB_B_aa(r)+sum_element_B_aa;
    m_B_aa=sum_element_B_aa/64;
end
for r=1:64
    diff_element_B_aa=(aveEn_dB_B_aa(r)-m_B_aa)^2+diff_element_B_aa;
    v_B_aa=diff_element_B_aa/64;
    sv_B_aa=v_B_aa^(0.5);
end
for r=1:64
    gPDF_B_aa(r)=(1/sqrt(2*pi*v_B_aa))*exp(-(aveEn_dB_B_aa(r)-m_B_aa)^2/2/v_B_aa);
end

gPDF_s(:,1)=(gPDF_A_s(1:80))';
gPDF_s(:,2)=(gPDF_B_s(1:80))';
gPDF_aa(:,1)=(gPDF_A_aa(1:64))';
gPDF_aa(:,2)=(gPDF_B_aa(1:64))';

% Making a table with the table function
item={'Region A of aa';'Region B of aa';'Region A of s';'Region B of s'};
mean={m_A_aa;m_B_aa;m_A_s;m_B_s};
variance={v_A_aa;v_B_aa;v_A_s;v_B_s};
table_para=table(item,mean,variance);

save('table with modelling of the energy values.mat','gPDF_aa','gPDF_s');
load ('table with modelling of the energy values.mat');

% figure;
% x_A_s=m_A_s-3*sv_A_s:0.01:m_A_s+3*sv_A_s;
% p_A_s=GPDF(x_A_s,m_A_s,v_A_s);
% plot(x_A_s,p_A_s);
% grid; xlabel('aveEn_dB_A_s'); ylabel('PGDF(aveEn_dB_A_s)');
% table(1,:)=p_A_s;

function [plot]=GPDF(x,mu,var)
z=size(x);
plot=ones(z)*(1/sqrt(2*pi*var));
me=ones(z)*mu;
v=ones(z)*var;
plot=plot.*(exp(-((x-me).*(x-me))./(2*v)));
end