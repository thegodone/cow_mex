clear all;
load X;
P = X(1,:);
T = X(2,:);

m=100;
t=20;
parfor i=1:100

tic
[warping1,Xw_cow1,diagnos1] = cowc_use_mex(P,T,m,t);
if mod(i,10)==0
    toc
end
end
