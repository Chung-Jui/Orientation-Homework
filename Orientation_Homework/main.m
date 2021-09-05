clc
clear all

%Parameter setting
F = [0;0;0;-1*10^7;0;0;0;-1*10^7];
sigma_y = 250*10^6;
E = 200*10^9;
L = 9.14;
density = 7860;

%object function and constrain
f = @(X) obj_fun(X,L,density);
g1 = @(X) abs(stress(X,F,E,L))-sigma_y;
g2 = @(X) node2_dis(X,F,E,L)-0.02;


x0 = [0.4 0.4];%initial design variables
lb = [0.001;0.001];%low bound
ub = [0.5;0.5];%upper bound
nonlcon = @(X) g_and_h(X,g1,g2);%constrain
A = []; 
b = [];
Aeq = [];
beq = [];
history = [0 x0 nonlcon(x0) f(x0)];
i=1;
[r_final fval history] = myproblem(f,x0,A,b,Aeq,beq,lb,ub,nonlcon,history,i);
[Q_final,sigma_final,R_final] = FEA(r_final,F,E,L)


%===================================plot================================================
fun = @(r1,r2) plot_fun(r1,r2,L,density);
fcontour(fun)
xlim([0.001 0.5])
ylim([0.001 0.5])
hold on
R=[];
for r1=-0.1:0.013:0.6
    for r2=-0.1:0.013:0.6
        r=[r1 r2];
        if g2(r)>0
            R = [R;r1 r2];
        end
    end
end
plot(R(:,1),R(:,2),'bx','MarkerSize',5)

R2=[];
for r1=-0.1:0.012:0.6
    for r2=-0.1:0.012:0.6
        r=[r1 r2];
        if max(g1(r)>0)
            R2 = [R2;r1 r2];
        end
    end
end
plot(R2(:,1),R2(:,2),'c*','MarkerSize',4)


X1 = history(:,2);X2 = history(:,3);

scatter3(X1,X2,zeros(length(X1),1))
plot(X1(end),X2(end),'ro','MarkerSize',10,'MarkerFaceColor',[1,0,0])
plot3(X1,X2,zeros(length(X1),1))

title('Matlab mathod(fmincon)','fontsize',15)
xlabel('r1 (m)','fontsize',15);
ylabel('r2 (m)  ','fontsize',15);
legend('Weight','£Gs_2 > 0.02(m)','¡U£m_i¡U> £m_y','optimal process','optimal point','Position',[0.75 0.85 0.12 0.05])
%=======================================================================================================

function [c,ceq] = g_and_h(X,g1,g2)
c = [g1(X) g2(X)];
ceq = [];
end








