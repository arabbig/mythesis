[JYVIXIndex, Date] = xlsread('implied.xlsx','W3:X1862'); 
JYVIXtab = table(Date,JYVIXIndex,'VariableNames',{'Date','Value'});
