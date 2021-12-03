N = 25;
s = tf('s');
gp = zpk(1/(s+1)); for i=2:N, gp = gp/(s+i); end, gp
pe = eig(gp); pp = eig(tf(gp));
plot(real(pe),imag(pe),'gx',real(pp),imag(pp),'r*')
title('Poles of 1/((s+1)(s+2)\cdot\cdot\cdot(s+25))')
ylabel('Imaginary Axis (seconds^{-1})')
xlabel('Real Axis (seconds^{-1})')
grid
legend('True poles','Perturbed poles')

