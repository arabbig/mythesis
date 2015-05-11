[pks locs] = findpeaks(coefs_eta, (0:LEN-1)');

figure;
subplot(2,1,1);
plot(index);
hold on
for i = 1:length(locs)
    d = locs(i) + 1;
    cumret = 1;
    while index(d)/index(d-1) < 1.03 && d < LEN
        if d - locs(i) > 10
            break
        end
        cumret = cumret * index(d)/index(d-1);
        d = d + 1;
    end
    if pks(i) > 0.2
        fprintf('eta=%6.3f days=%2.d drawdown=%8.3f%% \n', pks(i), d - locs(i) - 1,(cumret-1)*100);
        if cumret-1 > 0
            col = 'r';
        else
            col = 'g';
        end
        H=area(locs(i):d-1,index(locs(i):d-1));
        set(H,'FaceColor',col);
    end
end

subplot(2,1,2)
plot((0:LEN-1)',coefs_eta,locs,pks,'b*');