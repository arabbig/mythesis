require(tikzDevice)
seed(1234)
tikz('images/swarEx.tex',
     width=3.5,height=3.5)
x1 = arima.sim(n = 60, list(ar = c(-0.7)), sd = sqrt(0.0001))
x2 =  arima.sim(n = 60, list(ar = c(-0.1)), sd = sqrt(0.0001))
x <- c(x1,x2)
par(mar=c(3.5, 4.1, 2.1, 2.1), mgp=c(2.5, 1, 0), las=0)
plot(x, ylim = c(-0.04,0.065), xlab = 't', type='l', ylab ='y(t)')
text(c(30,90),c(0.06,0.04),labels=c('$y_t = -0.5 y_{t-1} + \\epsilon$','$y_t = -0.1 y_{t-1} + \\epsilon$'))
abline(v = 60, col='black', lty = 'dashed');
dev.off()

tikz('images/swarEx2.tex',  width=3.5,height=3.5)
s = rep(c(1,2),c(60,60));
plot(s,type='l',yaxt = 'n',ylab = 's(t)', xlab = 't');
axis(side = 2, at = c(1,2))
dev.off();