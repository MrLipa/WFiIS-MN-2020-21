a=VarName1;
b=VarName2;

t = 0:0.01:(1);
x = cos(t);

plot(a,b,"+",a,b,t,x);
title('Wychylenie x(t) oscylatora');
legend('x(t), dt=0.1','x(t)','cos(t)');
xlabel("t");
ylabel("x(t)");
grid();