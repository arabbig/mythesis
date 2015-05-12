function caldrawdown(fts_window, coefs_zeta, minzeta, fid)

index = fts2mat(fts_window.series1);
LEN = length(index);
[pks, locs] = findpeaks(coefs_zeta, (1:LEN)');

figure;
subplot(2,1,1);
plot(index);
hold on

cnt = 0;
correct = 0;
for i = 1:length(locs)
    d = locs(i) + 1;
    cumret = 1;
    while index(d)/index(d-1) < 1.029 && d < LEN
        if d - locs(i) > 10
            break
        end
        cumret = cumret * index(d)/index(d-1);
        d = d + 1;
    end
    if pks(i) > minzeta
        cnt = cnt + 1;
        fprintf(fid,'zeta=%6.3f days=%2.d drawdown=%8.3f%% \n', pks(i), d - locs(i) - 1,(cumret-1)*100);
        if cumret-1 > 0
            col = 'r';
        else
            correct = correct + 1;
            col = 'g';
        end
        H=area(locs(i):d-1,index(locs(i):d-1));
        set(H,'FaceColor',col);
    end
end

fprintf(fid,'Num of events=%d\nNegative Drawdown=%d\nP(Drawdown=True)=%f \n', cnt, correct, correct/cnt);

subplot(2,1,2)
%plot((0:LEN-1)',coefs_zeta,locs,pks,'b*');
[ax,h1,h2] = plotyy((1:LEN),index,(1:LEN),coefs_zeta);
set(ax(2),'YLim',[0 0.5]);
end