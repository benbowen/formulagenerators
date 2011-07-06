clear all
close all
clc
%%
mex HR2_matlab_str.cpp
%%
tic
a=HR2_matlab_str('-m' ,'544.142307' ,'-t' ,'2.7' ,'-C' ,'0-200' ,'-H' ,'0-400' ,'-N' ,'0-10' ,'-O' ,'0-10' ,'-P' ,'0-6', '-S', '0-4', '-F', '0-0', '-L' ,'0-0' ,'-B' ,'0-0' ,'-I','0-0');
formula=[]
toc
for i = 1:size(a,1);
    formula{i}=FormulaToStructure(a{i});
end
toc

%%
mex HR2_bruteForce_matlab_str.cpp
tic
a=HR2_bruteForce_matlab_str('-m' ,'544.142307' ,'-t' ,'2.7' ,'-C' ,'0-200' ,'-H' ,'0-400' ,'-N' ,'0-10' ,'-O' ,'0-10' ,'-P' ,'0-6', '-S', '0-4', '-F', '0-0', '-L' ,'0-0' ,'-B' ,'0-0' ,'-I','0-0');
toc
size(a)