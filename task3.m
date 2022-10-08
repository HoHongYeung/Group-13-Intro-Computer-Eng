clc
clear all
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
        timeSegStart_ms=str2num(aasc1)*10^(-6)+(str2num(aasc2)*10^(-6)-str2num(aasc1)*10^(-6))/2-12.5;
        timeSegEnd_ms=str2num(aasc2)*10^(-6)-((str2num(aasc2)*10^(-6)-str2num(aasc1)*10^(-6))/2-12.5);
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
        timeSegStart_ms=str2num(aasc1)*10^(-6)+(str2num(aasc2)*10^(-6)-str2num(aasc1)*10^(-6))/2-12.5;
        timeSegEnd_ms=str2num(aasc2)*10^(-6)-((str2num(aasc2)*10^(-6)-str2num(aasc1)*10^(-6))/2-12.5);
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
% segOrig_phAA(1,:);

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
