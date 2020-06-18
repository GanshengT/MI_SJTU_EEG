%% 使用libsvm 运动想象分类
clc;

disp('#######  Training The SVM Classsifier ##########')
load dataCSP_1.mat 
load test_lab.mat

%% X为训练集， T测试集， Y训练标签

%% SVM网络训练
model = libsvmtrain(Y, X, '-c 2 -g 1'); %fitcsvm

%% SVM网络预测
[predict_label, accuracy,decision_values] = libsvmpredict(y_test, T, model);
