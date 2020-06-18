% Extract Common Spatial Pattern (CSP) Feature
close all; clear; clc;

% bb=load('F:\Matlabcx\CPS_zzl_0921_train');	
bb=load('train_zzl');	
% cc=load('F:\Matlabcx\CPS_zzl_0921_test');	
cc=load('test_zzl');	
trian_label=load('train_label.txt');	
y_test=load('test_label.txt');	
%y_test=trian_label(51:100,:);	
EEG_tr=bb(2:33,:);
EEG_tr=EEG_tr';
EEG_te=bb(2:33,:);
EEG_te=EEG_te';
EEG_train=zeros(2001,32,50);
EEG_test=zeros(2001,32,50);
%for i=1:50
%    EEG_train(:,:,i)=EEG_tr(1501+(i-1)*4000:3501+(i-1)*4000,:);
%end
%for i=1:50
%    EEG_test(:,:,i)=EEG_te(1501+(i+50-1)*4000:3501+(50+i-1)*4000,:);
%end

for i=1:100
    EEG_train(:,:,i)=EEG_tr(1501+(i-1)*4000:3501+(i-1)*4000,:);
    EEG_test(:,:,i)=EEG_te(1501+(i-1)*4000:3501+(i-1)*4000,:);
end

EEGSignals.x=EEG_train;
EEGSignals.y=trian_label;
Y=trian_label;

classLabels = unique(EEGSignals.y); 
CSPMatrix = learnCSP(EEGSignals,classLabels);
nbFilterPairs = 1;

X = extractCSP(EEGSignals, CSPMatrix, nbFilterPairs);  
% EEGSignals.x=EEG_train(:,:,51:100);
EEGSignals.x=EEG_test;
%EEGSignals.y=y_test;%这一行是后加的代码，是否正确还不确定
%CSPMatrix2 = learnCSP(EEGSignals,classLabels);%这一行是后加的代码，是否正确还不确定
T = extractCSP(EEGSignals, CSPMatrix, nbFilterPairs); 

% [inputn,inputps]=mapminmax(X);
% X=inputn;
% [inputn_T,inputps_T]=mapminmax(T);
% T=inputn_T;

save dataCSP_1.mat X Y T
save test_lab.mat y_test

color_L = [0 102 255] ./ 255;
color_R = [255, 0, 102] ./ 255;
pos = find(Y==1);
plot(X(pos,1),X(pos,2),'x','Color',color_L,'LineWidth',2);
hold on
pos = find(Y==2);
plot(X(pos,1),X(pos,2),'o','Color',color_R,'LineWidth',2);
legend('Left Hand','Right Hand')
xlabel('C3','fontweight','bold')
ylabel('C4','fontweight','bold')