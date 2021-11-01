% Inflation target plotter (credibility example)

for t=1:T_sim-1

        pi(t) = X_stack(1,t); y(t) = X_stack(2,t); int(t) = X_stack(3,t); 
        
end 

set(0,'DefaultLineLineWidth',1.5)
hold on, 
subplot(2,3,4), hold on, plot(Periods, 4*100*pi, 'g'), hold on, title('Inflation'), xlabel('Periods'), ylabel('% (annualized)'), axis([-inf,inf,-inf,inf]) %axis([1,16,2,6])
subplot(2,3,5), hold on, plot(Periods, 100*y,'g'), title('Output'), xlabel('Periods'), ylabel('%'), axis([-inf,inf,-inf,inf]), %axis([1,16,-0.2,0.6])
subplot(2,3,6), hold on, plot(Periods, 4*100*int,'g'), hold on, title('Nominal interest rate'), axis([-inf,inf,-inf,inf])  %axis([1,16,5,9.01]),
xlabel('Periods'), ylabel('% (annualized)'), hold on