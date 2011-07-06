function c = recursiveloop(bnds)
c=1;
for i = 1:size(bnds,2)
    c=otherloop(c,bnds(1,i),bnds(2,i))
end

function c = otherloop(c,lb,ub)
    for i = lb:ub
        c=c+1;
    end


