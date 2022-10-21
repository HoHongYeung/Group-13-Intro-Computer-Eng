clear;
disp('Welcome to a text-based menu-driven program for calculation of PI');
disp('Press a -- Perform FIR filtering');
disp('Press b -- Extract signal segments');
disp('Press c -- Calculate energy for specified frequency regions');
disp('Press d -- Modelling of energy values using Gaussian PDF');
disp('Press e -- Exit the program');
strResponse = input('Please make your choice: ', 's');
while (strcmp(strResponse,'e')~=1)
    switch (strResponse)
        case 'a'
            fid = fopen('C:\Users\SHLT\Desktop\listData.txt','rt');
            textLine1 = fscanf(fid,'%s',1);
            cishu=0;
            while (strcmp(textLine1,'.')~=1)
            %   textLine1 = fscanf(fid,'%s',1);
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
  if(cishu==0)
    figure;
    plot([1900:2300],inpSigWav(1900:2300),'g',[1900:2300],outSigWav(1900:2300),'b'); 
    grid; xlabel('Sample index'); ylabel('Amplitude');
    legend('Input signal', 'Output signal')
%     bb(1)=upper(bb(1));
%     bb(2)=upper(bb(2));
  end
  fileNameOut = strcat('C:\Users\SHLT\Desktop\wavFilt','\',aa,'\',bb);
  audiowrite(fileNameOut, outSigWav, Fs);
  textLine1 = fscanf(fid,'%s',1);

  cishu=cishu+1;
end
fupintexing();

        case 'b'
           k=0;
numaa=0;
nums=0;
fid = fopen('C:\Users\SHLT\Desktop\listData.txt','rt');
textLine1 = fscanf(fid,'%s',1);
aa=textLine1(1:5);
bb=textLine1(7:end);
bb1=strcat(bb,'.lab');
fileName=strcat('C:\Users\SHLT\Downloads\dataTIMIT_labAssign2022_usedToStud\labels','\',aa,'\',bb1);
fidd=fopen(fileName);
%读lab
while (strcmp(textLine1,'.')~=1)
    
    aasc1=fscanf(fidd,'%s',1);
    aasc2=fscanf(fidd,'%s',1);
    aascx=fscanf(fidd,'%s',1);
    if(strcmp(aascx,'aa')==1)
        numaa=numaa+1;
        timeSegStart_ms=str2num(aasc1)*10^(-4)+(str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5;
        timeSegEnd_ms=str2num(aasc2)*10^(-4)-((str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5);
%         if(timeSegStart_ms<0)
%             timeSegStart_ms=0;
%         end
        bb(1)=lower(bb(1));
        bb(2)=lower(bb(2));
        bb=bb(1:end);
        bb2=strcat(bb,'.wav');
        fileName=strcat('C:\Users\SHLT\Desktop\wavOrig','\',aa,'\',bb2);
        [inpSigWav,Fs] = audioread(fileName);
        diankaishi=timeSegStart_ms*10^(-3)/(1/8000);
        if(diankaishi<0)
            diankaishi=1;
        end
        diankaishi=ceil(diankaishi);
        dianjieshu=timeSegEnd_ms*10^(-3)/(1/8000);
        dianjieshu=ceil(dianjieshu);
        changdu=length(inpSigWav(diankaishi:dianjieshu));
        segOrig_phAA(numaa,:)=changdu;
        if(numaa==1)
            segOrig_phAA(1,:)
        end
    end
    if(strcmp(aascx,'s')==1)
        nums=nums+1;
        timeSegStart_ms=str2num(aasc1)*10^(-4)+(str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5;
        timeSegEnd_ms=str2num(aasc2)*10^(-4)-((str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5);
%         if(timeSegStart_ms<0)
%             timeSegStart_ms=0;
%         end
        bb(1)=lower(bb(1));
        bb(2)=lower(bb(2));
        bb=bb(1:end);
        bb2=strcat(bb,'.wav');
        fileName=strcat('C:\Users\SHLT\Desktop\wavOrig','\',aa,'\',bb2);
        [inpSigWav,Fs] = audioread(fileName);
        diankaishi=timeSegStart_ms*10^(-3)/(1/8000);
        if(diankaishi<0)
            diankaishi=1;
        end
        diankaishi=ceil(diankaishi);
        dianjieshu=timeSegEnd_ms*10^(-3)/(1/8000);
        dianjieshu=ceil(dianjieshu);
        changdu=length(inpSigWav(diankaishi:dianjieshu));
        segOrig_phS(nums,:)=changdu;
    end
   if(strcmp(aasc1,'')==1)
      textLine1 = fscanf(fid,'%s',1);
      if(strcmp(textLine1,'.')==0)
        aa=textLine1(1:5);
        bb=textLine1(7:end);
        bb1=strcat(bb,'.lab');
        fileName=strcat('C:\Users\SHLT\Downloads\dataTIMIT_labAssign2022_usedToStud\labels','\',aa,'\',bb1);
        fidd=fopen(fileName);
      end
   end
end
save segAllData segOrig_phAA segOrig_phS;
% segOrig_phAA(1,:);
        case 'c'
            yinsu={'aa';'s'};
k=0;
numaa=0;
nums=0;
fid = fopen('C:\Users\SHLT\Desktop\listData.txt','rt');
textLine1 = fscanf(fid,'%s',1);
aa=textLine1(1:5);
bb=textLine1(7:end);
bb1=strcat(bb,'.lab');
fileName=strcat('C:\Users\SHLT\Downloads\dataTIMIT_labAssign2022_usedToStud\labels','\',aa,'\',bb1);
fidd=fopen(fileName);
%读lab
while (strcmp(textLine1,'.')~=1)
    
    aasc1=fscanf(fidd,'%s',1);
    aasc2=fscanf(fidd,'%s',1);
    aascx=fscanf(fidd,'%s',1);
    if(strcmp(aascx,'aa')==1)
        numaa=numaa+1;
        timeSegStart_ms=str2num(aasc1)*10^(-4)+(str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5;
        timeSegEnd_ms=str2num(aasc2)*10^(-4)-((str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5);
%         if(timeSegStart_ms<0)
%             timeSegStart_ms=0;
%         end
        bb(1)=lower(bb(1));
        bb(2)=lower(bb(2));
        bb=bb(1:end);
        bb2=strcat(bb,'.wav');
        fileName=strcat('C:\Users\SHLT\Desktop\wavOrig','\',aa,'\',bb2);
        [inpSigWav,Fs] = audioread(fileName);
        diankaishi=timeSegStart_ms*10^(-3)/(1/8000);
        if(diankaishi<0)
            diankaishi=1;
        end
        diankaishi=ceil(diankaishi);
        dianjieshu=timeSegEnd_ms*10^(-3)/(1/8000);
        dianjieshu=ceil(dianjieshu);
        changdu=length(inpSigWav(diankaishi:dianjieshu));


        fFn=8000/changdu;
        kai=400/fFn;
        kai=ceil(kai)+1;
        jie=1600/fFn;
        jie=ceil(jie)+1;
        zuizhong=jsdft(inpSigWav(diankaishi:dianjieshu),changdu);
        yy=jsnl(zuizhong,kai,jie);
        enOrigRegAB_phAA(numaa,1)=yy;

        kai=2400/fFn;
        kai=ceil(kai)+1;
        jie=4000/fFn;
        jie=ceil(jie)+1;
        zuizhong=jsdft(inpSigWav(diankaishi:dianjieshu),changdu);
        yy=jsnl(zuizhong,kai,jie);
        enOrigRegAB_phAA(numaa,2)=yy;



%       segOrig_phAA(numaa,:)=changdu;
%         if(numaa==1)
%             segOrig_phAA(1,:)
%         end
    end
    if(strcmp(aascx,'s')==1)
        nums=nums+1;
        timeSegStart_ms=str2num(aasc1)*10^(-4)+(str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5;
        timeSegEnd_ms=str2num(aasc2)*10^(-4)-((str2num(aasc2)*10^(-4)-str2num(aasc1)*10^(-4))/2-12.5);
%         if(timeSegStart_ms<0)
%             timeSegStart_ms=0;
%         end
        bb(1)=lower(bb(1));
        bb(2)=lower(bb(2));
        bb=bb(1:end);
        bb2=strcat(bb,'.wav');
        fileName=strcat('C:\Users\SHLT\Desktop\wavOrig','\',aa,'\',bb2);
        [inpSigWav,Fs] = audioread(fileName);
        diankaishi=timeSegStart_ms*10^(-3)/(1/8000);
        if(diankaishi<0)
            diankaishi=1;
        end
        diankaishi=ceil(diankaishi);
        dianjieshu=timeSegEnd_ms*10^(-3)/(1/8000);
        dianjieshu=ceil(dianjieshu);
        changdu=length(inpSigWav(diankaishi:dianjieshu));
        
        fFn=8000/changdu;
        kai=400/fFn;
        kai=ceil(kai)+1;
        jie=1600/fFn;
        jie=ceil(jie)+1;
        zuizhong=jsdft(inpSigWav(diankaishi:dianjieshu),changdu);
        yy=jsnl(zuizhong,kai,jie);
        enOrigRegAB_phS(nums,1)=yy;

        kai=2400/fFn;
        kai=ceil(kai)+1;
        jie=4000/fFn;
        jie=ceil(jie)+1;
        zuizhong=jsdft(inpSigWav(diankaishi:dianjieshu),changdu);
        yy=jsnl(zuizhong,kai,jie);
        enOrigRegAB_phS(nums,2)=yy;

    end
   if(strcmp(aasc1,'')==1)
      textLine1 = fscanf(fid,'%s',1);
      if(strcmp(textLine1,'.')==0)
        aa=textLine1(1:5);
        bb=textLine1(7:end);
        bb1=strcat(bb,'.lab');
        fileName=strcat('C:\Users\SHLT\Downloads\dataTIMIT_labAssign2022_usedToStud\labels','\',aa,'\',bb1);
        fidd=fopen(fileName);
      end
   end
end
energyA={enOrigRegAB_phAA(1,1);enOrigRegAB_phS(1,1)};
energyB={enOrigRegAB_phAA(1,2);enOrigRegAB_phS(1,2)};
T=table(yinsu,energyA,energyB)
figure(1);
T1=histogram(enOrigRegAB_phAA(:,1));
figure(2);
T2=histogram(enOrigRegAB_phAA(:,2));
figure(3);
T3=histogram(enOrigRegAB_phS(:,1));
figure(4);
T4=histogram(enOrigRegAB_phS(:,2));
% segOrig_phAA(1,:);

        case 'd'
            fid = fopen('C:\Users\SHLT\Desktop\listData.txt','rt');
textLine = fscanf(fid,'%s',1);
ind_s = 0;
ind_aa = 0;
num_s=0;
num_aa=0;
Orig_s=[];
Orig_aa=[];%Store full audio

while (strcmp(textLine,'.')~=1)
    filename=strcat('C:\Users\SHLT\Downloads\dataTIMIT_labAssign2022_usedToStud\labels\',textLine,'.lab');
    file=fopen(filename,'rt');
    textCol = textscan(file,'%f %f %s');%Scan the 3 columns
    audioname=strcat('C:\Users\SHLT\Desktop\wavOrig\',textLine,'.wav');
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
    xxx(j,1)=0;
    xxx(j,2)=0;%for frequency region B
    %intercept the Orig_s according to the time index, and store 
    %the result in segOrig_phS
    segOrig_phS(j,1:201)=(Orig_s{j}(indexSegStart_s(j):indexSegEnd_s(j),1))';
    segOrig_phS_ft(j,1:201)=abs(fft(segOrig_phS(j,1:201)));
    %the |X(k)|^2 of the number of columns 11 to 41 of 
    %each row is added and stored into the first column of value1
    for m=ind1_A:ind2_A
        xxx(j,1)=abs(segOrig_phS_ft(j,m))^2+xxx(j,1);
    end
    aveEn_dB_A_s(j,1)=10*log10(xxx(j,1)/(ind2_A-ind1_A+1));
    enOrigRegAB_phS(j,1)=aveEn_dB_A_s(j,1);
    for n=ind1_B:ind2_B
        xxx(j,2)=abs(segOrig_phS_ft(j,n))^2+xxx(j,2);
    end
    aveEn_dB_B_s(j,1)=10*log10(xxx(j,2)/(ind2_B-ind1_B+1));
    enOrigRegAB_phS(j,2)=aveEn_dB_B_s(j,1);
end

for k=1:ind_aa
    yyy(k,1)=0;
    yyy(k,2)=0;
    segOrig_phAA(k,1:201)=(Orig_aa{k}(indexSegStart_aa(k):indexSegEnd_aa(k),1))';
    segOrig_phAA(k,1:201)=(Orig_aa{k}(indexSegStart_aa(k):indexSegEnd_aa(k),1))';
    segOrig_phAA_ft(k,1:201)=abs(fft(segOrig_phAA(k,1:201)));
    for o=ind1_A:ind2_A
        yyy(k,1)=abs(segOrig_phAA_ft(k,o))^2+yyy(k,1);
    end
    aveEn_dB_A_aa(k,1)=10*log10(yyy(k,1)/(ind2_A-ind1_A+1));
    enOrigRegAB_phAA(k,1)=aveEn_dB_A_aa(k,1);
    for p=ind1_B:ind2_B
        yyy(k,2)=abs(segOrig_phAA_ft(k,p))^2+yyy(k,2);
    end
    aveEn_dB_B_aa(k,1)=10*log10(yyy(k,2)/(ind2_B-ind1_B+1));
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
meannn={m_A_aa;m_B_aa;m_A_s;m_B_s};
varianceee={v_A_aa;v_B_aa;v_A_s;v_B_s};
table_para=table(item,meannn,varianceee);

save('table with modelling of the energy values.mat','gPDF_aa','gPDF_s');
load ('table with modelling of the energy values.mat');

% figure;
% x_A_s=m_A_s-3*sv_A_s:0.01:m_A_s+3*sv_A_s;
% p_A_s=GPDF(x_A_s,m_A_s,v_A_s);
% plot(x_A_s,p_A_s);
% grid; xlabel('aveEn_dB_A_s'); ylabel('PGDF(aveEn_dB_A_s)');
% table(1,:)=p_A_s;

        otherwise
           disp('incorrect choice - please enter again'); 
    end
    strResponse = input('Please make your choice: ', 's');
end

function[y]=jsdft(h,l)
%l是总点数
  y=fft(h,l);
end

function[yy]=jsnl(k,i1,i2)
%k为计算出来的DFT
zz=0;
for i=i1:i2
    zz=abs(k(i))^2+zz;
end
 yy=10*log10(1/(i2-i1+1)*zz);
end

function [plot]=GPDF(x,mu,var)
z=size(x);
plot=ones(z)*(1/sqrt(2*pi*var));
me=ones(z)*mu;
v=ones(z)*var;
plot=plot.*(exp(-((x-me).*(x-me))./(2*v)));
end
