%% ʹ��libsvm �˶��������
clc;

disp('#######  Training The SVM Classsifier ##########')
load dataCSP_1.mat 
load test_lab.mat

%% XΪѵ������ T���Լ��� Yѵ����ǩ

%% SVM����ѵ��
model = libsvmtrain(Y, X, '-c 2 -g 1'); %fitcsvm

%% SVM����Ԥ��
[predict_label, accuracy,decision_values] = libsvmpredict(y_test, T, model);
