%% Analysis Script
clear all

%Specify variables FIRST
FOLDER='E.coli30r30s/';
OVERLAP_line=[20,40,60,80,100,120,200,500,1000];

DATASET_size=30;
READ_length=5000;
THRESHOLD=0.9;
ERR=1;

%Generate results tables
Results=cell(length(OVERLAP_line)+1,7);
Results{1,1}='Overlaps length';Results{1,2}='Accuracy';Results{1,3}='Recall/Sensitivity';Results{1,4}='Precision/PPV';
Results{1,5}='Correct joins';Results{1,6}='Incorrect joins';Results{1,7}='Missing';
%Generate table Final results
ResultsX=cell(length(OVERLAP_line)+1,7);ResultsX{1,1}='Overlaps length';ResultsX{1,2}='Accuracy';ResultsX{1,3}='Recall/Sensitivity';ResultsX{1,4}='Precision/PPV';
ResultsX{1,5}='Correct joins %';ResultsX{1,6}='Incorrect joins %';ResultsX{1,7}='Missing';
for i=1:length(OVERLAP_line)
    Results{i+1,1}=OVERLAP_line(i);
    ResultsX{i+1,1}=OVERLAP_line(i);
end

%% Running and saving results
ACC=NaN(1,DATASET_size);SSV=NaN(1,DATASET_size);PPV=NaN(1,DATASET_size);JNS=[];DJNS=[];MSNG=[];

for j=1:length(OVERLAP_line)
    for i=1:DATASET_size
        filename=([FOLDER 'reads ' num2str(ERR) '% errors ' num2str(READ_length) '-' num2str(OVERLAP_line(j)) ' set' num2str(i) '.fasta']);
%         disp(filename)        
        [header,fasta]=fastaread(filename);
        reference=[1:length(fasta)];
        
        result=GS_overlap_detection(fasta,THRESHOLD);
        
        disp(['OK ' num2str(READ_length) '-' num2str(OVERLAP_line(j)) ' set' num2str(i)])
        [accuracy,sensitivity,positive_predictivity,joins,disjoins,missing]=GS_overlap_detection_results(reference,result);
        ACC(i)=accuracy;
        SSV(i)=sensitivity;
        PPV(i)=positive_predictivity;
        JNS(i)=joins;
        DJNS(i)=disjoins;
        MSNG(i)=missing;
        
    end
    Results{j+1,2}=ACC;
    Results{j+1,3}=SSV;
    Results{j+1,4}=PPV;
    Results{j+1,5}=JNS;
    Results{j+1,6}=DJNS;
    Results{j+1,7}=MSNG;
    ResultsX{j+1,2}=nanmean(ACC);
    ResultsX{j+1,3}=nanmean(SSV);
    ResultsX{j+1,4}=nanmean(PPV);
    
    ResultsX{j+1,5}=mode(JNS);
    ResultsX{j+1,6}=mode(DJNS);
    ResultsX{j+1,7}=mode(MSNG);
end
format shortg
c=clock;
save(['Results ' num2str(READ_length) 'bp ' num2str(ERR) '% errors ' num2str(THRESHOLD) 'threshold ' num2str(c(4)) '.' num2str(c(5)) '.mat'],'Results','ResultsX')
