SVMModel = fitcsvm(X,Y,'kernelFunction','gaussian')
[label,score] = predict(SVMModel,T)
accuracy = 1-mean(label - y_test)