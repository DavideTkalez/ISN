function [acc,spec,sen] = calcolo_metriche(CM)

TP = CM(1,1);
TN = CM(2,2);
FP = CM(1,2);
FN = CM(2,1);

acc = (TP+TN)/sum(CM,"all");

spec = TN/(TN + FP);
sen = TP/(TP+FN);


end