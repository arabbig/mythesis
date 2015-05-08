n = 0:319;
x = cos(pi/4*n)+0.5*sin(pi/2*n)+randn(size(n));

[pxx1,w1] = periodogram(x);
plot(w1/pi,pxx1)
legend('pxx1')
xlabel('\omega / \pi')