clear;
disp('Welcome to a text-based menu-driven program for analysis of a given set of audio data');
disp('Press a -- Perform FIR filtering');
disp('Press b -- Extraction of signal segments');
disp('Press c -- Calculate the average energy for specified frequency regions');
disp('Press d -- Modelling of the energy values using Gaussian PDF');
disp('Press e -- Exit the program');
strResponse = input('Please make your choice: ', 's');
while (strcmp(strResponse,'e')~=1)
    switch (strResponse)
        case 'a'
            clear;
            fid = fopen('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\listData.txt','rt');
            ind=0;            
            % Read the column(column 1)
            textCol1 = fscanf(fid,'%s',1);
            while (strcmp(textCol1,'.')~=1)%While the content of column 1 is not '.'
                ind = ind+1;
                save_path='C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavOrig\';
                content='textLine1(1:6).textLine1(7:8).textLine1(9:end)';
                cell_str=strsplit(content,'.');
                ucl=textCol1(7:8);
                if ucl>='A' & ucl<='Z'
                    lcl=lower(ucl);
                end        
                sep_name=cell_str{1,1};
                newname=strcat(save_path,textCol1(1:6),lcl,textCol1(9:end),'.wav');
                
                % Store the data into a structure array
                storeDataStruct(ind).name = newname;              
                
                [inpSigWav,Fs] = audioread(newname);%  Load the input signal
                nSamples = length(inpSigWav);
                
                %  Define the FIR filter coefficients
                firCoef = [-0.08,0.24,0.4,0.4,0.16,-0.24,0.08];
                nCoef = length(firCoef);
                
                %  Add zero samples to the beginning and to the end of the input signal to take account of the filter length
                inpSigWavExt = [zeros(6,1); inpSigWav; zeros(6,1)];
                nSamplesNew = length(inpSigWavExt);
                
                %  Allocate memory for the output signal
                outSigWav = zeros(nSamplesNew,1);
                
                %  Perform filtering
                for i=nCoef:nSamplesNew,
                    outSigWav(i) = firCoef(1)*inpSigWavExt(i)+firCoef(2)*inpSigWavExt(i-1)+firCoef(3)*inpSigWavExt(i-2)+firCoef(4)*inpSigWavExt(i-3)+firCoef(5)*inpSigWavExt(i-4)+firCoef(6)*inpSigWavExt(i-5)+firCoef(7)*inpSigWavExt(i-6);
                end
                outSigWav=outSigWav(7:end-6);%Should it be?
                
                
                %  Store the output signal into .wav file
                fileNameOut = strcat('C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavFilt\',textCol1(1:6),lcl,textCol1(9:end),'.wav');
                audiowrite(fileNameOut, outSigWav, Fs);
                if (strcmp(textCol1,'MCPM0\SA1')==1)%test
                    x=outSigWav;
                end
                
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
            
            fileNameOrig = 'C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavOrig\MCPM0\sa1.wav';
            [inpSigWav,Fs] = audioread(fileNameOrig);%  Load the input signal
            
            fileNameFilt = 'C:\Users\28230\Desktop\dataTIMIT_labAssign2022_usedToStud\wavFilt\MCPM0\sa1.wav';
            [outSigWav,Fs] = audioread(fileNameFilt);%  Load the output signal
            
            figure;
            plot([1 : 401], inpSigWav(1900:2300),'r',[1 : 401],outSigWav(1900:2300),'g');
            % plot([1900 : 2300], x(1900:2300),'r',[1900 : 2300],outSigWav(1900:2300),'g');            
            grid; xlabel('Index'); ylabel('Amplitude');
            legend('Input signal', 'Output signal')
            
            %test
            % y=x(7:24377);
            % z=filter(firCoef,1,inpSigWav);
            
            fclose(fid);
        case 'b'
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
        case 'c'
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
            
            figure;
            enOrigRegA_phS=histogram(enOrigRegAB_phS(:,1));
            legend('the average energy values for the frequency region A for the phoneme ‘s’');
            figure;
            enOrigRegB_phS=histogram(enOrigRegAB_phS(1:end,2));
            legend('the average energy values for the frequency region B for the phoneme ‘s’');
            figure;
            enOrigRegA_phAA=histogram(enOrigRegAB_phAA(:,1));
            legend('the average energy values for the frequency region A for the phoneme ‘aa’');
            figure;
            enOrigRegB_phAA=histogram(enOrigRegAB_phAA(1:end,2));
            legend('the average energy values for the frequency region B for the phoneme ‘aa’');
            
            %Build the table
            
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
                for q=1:size(textCol1{1},1)
                    if strcmp(textCol1{3}(q),'s')==1
                        ind_s_sa1=ind_s_sa1+1;
                        if (strcmp(textLine1,'MDPK0\SA1')==1)
                            aa=ind_s_sa1;
                            aveEn_sa1_s=enOrigRegAB_phS(ind_s_sa1,1:end);
                            x1_s=indexSegStart_s(ind_s_sa1);
                            x2_s=indexSegEnd_s(ind_s_sa1);
                            break;
                        end
                    end
                end
                
                for r=1:size(textCol1{1},1)
                    if strcmp(textCol1{3}(r),'aa')==1
                        ind_aa_sa1=ind_aa_sa1+1;
                        if (strcmp(textLine1,'MDPK0\SA1')==1)
                            bb=ind_aa_sa1;
                            aveEn_sa1_aa=enOrigRegAB_phAA(ind_aa_sa1,1:end);
                            x1_aa=indexSegStart_aa(ind_aa_sa1);
                            x2_aa=indexSegEnd_aa(ind_aa_sa1);
                            break;
                        end
                    end
                end
                fclose(file1);
                textLine1 = fscanf(fid1,'%s',1);
            end
            
            save('table.mat','aveEn_sa1_aa','aveEn_sa1_s');
            load ('table.mat');
            
            fclose(fid1);
        case 'd'
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
            for q=1:80
                sum_element_A_s=aveEn_dB_A_s(q)+sum_element_A_s;
                m_A_s=sum_element_A_s/80;
            end
            for q=1:80
                diff_element_A_s=(aveEn_dB_A_s(q)-m_A_s)^2+diff_element_A_s;
                v_A_s=diff_element_A_s/80;
                sv_A_s=v_A_s^(0.5);
            end
            
            for q=1:80
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
        otherwise
            disp('incorrect choice - please enter again');
    end
    strResponse = input('Please make your choice: ', 's');
end
function [plot]=GPDF(x,mu,var)
            z=size(x);
            plot=ones(z)*(1/sqrt(2*pi*var));
            me=ones(z)*mu;
            v=ones(z)*var;
            plot=plot.*(exp(-((x-me).*(x-me))./(2*v)));
end
