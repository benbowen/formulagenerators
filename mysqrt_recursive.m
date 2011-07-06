function sq = mysqrt_recursive(x,g)
clc
if isempty(g)
    g = x/2;
end
t=0;
tic
while ~t
    [t,sq]=test(x,g);
end
toc


function [t,sq]=test(x,g)
if closeEnough(x/g, g)
    t=1;
    sq=g;
else
    [t,sq]=test(x,betterGuess(x,g));
end


function F=closeEnough(a,b)
F=abs(a-b)<1e-3;

function g=betterGuess(x,g)
g=(g+x/g)/2
