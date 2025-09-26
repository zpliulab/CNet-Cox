function [y0, pred_class, AUC, RefereceResult] = PredictSVM(train_data, train_y, Est_theta) 
%%
% Est_theta = theta;
% train_data = X_test;
% train_y = y_test;
[n1,p1] = size(train_data); 
X = [ones(n1,1), train_data]; 
y0 = train_y; 
theta = Est_theta;
%% 2 different data possible: for GENES and for CLASSES (here only for classes)
% separating line: x*w+b=f(x)
pp = X*theta;
% MATLAB has its own built-in function, plotroc, for plotting ROC curves with SVM classifiers, 
% but this function is not suitable for random forest classification models and predictions.
% AUC = plotroc(y0,pp);  
% AUC = AUC(y0,pp);
% pred_class = sign(pp);
pred_class = sign(pp);
figure;
AUC = plotROC(y0,pred_class);
%% 3. sensitivity and specificity  for CLASSES
tabClass = confusionmat(y0,pred_class);
figure;
tabPlot = printConMat(tabClass);
% stats = statsOfMeasure(tabClass);   % There are errors
[Result,RefereceResult] = confusion.getMatrix(y0,pred_class);
% disp(Result)
% disp(RefereceResult)
return