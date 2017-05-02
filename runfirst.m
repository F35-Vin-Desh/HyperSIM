%runfirst
val = xlsread('outputs.xlsx');
TIME = val(:,1);
TIME(length(TIME)) = [];
X = val(:,2);
X(length(X)) = [];
HEIGHT = val(:,3);
HEIGHT(length(HEIGHT))= [];
PITCH = val(:,4);
PITCH(length(PITCH)) = [];
%roll back from 180 to 0 (pitch up in): 360s
