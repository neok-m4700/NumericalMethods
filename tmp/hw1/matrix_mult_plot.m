load matrix_mult.dat;
x = matrix_mult(:,1);
m_loop = matrix_mult(:,2);
m_builtin = matrix_mult(:,3);
c_loop = matrix_mult(:,4);
figure
semilogy(x,m_loop, x,m_builtin, x,c_loop);
title('Semi log plot of running times for matrix multiplication vs size of side');
xlabel('side length')
ylabel('running time in seconds')
legend('matlab loop', 'matlab builtin', 'c loop')

