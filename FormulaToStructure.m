function formula = FormulaToStructure(str)
numpar=regexp(str,'[A-Z]');
for i = 1:length(numpar)
    if i<length(numpar)
        edges2(i)=numpar(i+1)-1;
    else
        edges2(i)=length(str);
    end
end

for e = 1:length(numpar)
    sym=regexp(str(numpar(e):edges2(e)),'[a-z_A-Z]+','match');
    num=regexp(str(numpar(e):edges2(e)),'\d+','match');
    if isempty(num)
        formula.(sym{1})=1;
    else
        formula.(sym{1})=str2double(num{1});
    end
end
