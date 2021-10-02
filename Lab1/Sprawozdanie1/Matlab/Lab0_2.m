a=VarName1;
b=e1;

plot(a,b,"+",a,b);
title('Błąd średniokwadratowy');
legend('o(q), q','o(q)');
xlabel("q");
ylabel("o(q)");
grid();