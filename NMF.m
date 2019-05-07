function [U1,U2,V1,V2,Vstar]=NMF(X1,X2,weight,K)

%行归一化矩阵X,viewNum即为算法中的对应视角
[cellNum1,viewNum]=size(X1);
sumX1=sum(X1,2);
D=diag(sumX1);
X1=D^-1*X1;

[cellNum2,viewNum]=size(X2);
sumX2=sum(X2,2);
D=diag(sumX2);
X2=D^-1*X2;

%初始化U，V，V*，为非负数
U1=rand(cellNum1,K);
V1=rand(viewNum,K);
U2=rand(cellNum2,K);
V2=rand(viewNum,K);
Vstar=rand(viewNum,K);

%记录o以便画图
iter=30;
plot_o1=zeros(1,iter);
plot_oo1=zeros(1,iter);
plot_o2=zeros(1,iter);
plot_oo2=zeros(1,iter);
plot_o=zeros(1,iter);
plot_oo=zeros(1,iter);

%M-NMF
Q1=[];
Q2=[];
v=zeros(1,K);
O=150;
O1=150;
O2=150;
O1_old=200;
O2_old=200;
O_old=200;
%while (abs(O_old-O)>1e-3)
for it=1:50
    for view=1:2
        %X1
        if view==1
            %while (abs(O1_old-O1)>1e-3)
            for it1=1:iter
                %Fixing Vstar&V1, update U1
                u_up=X1*V1;
                u_down=U1*(V1')*V1;
                for k=1:K
                    for i=1:cellNum1
                        U1(i,k)=U1(i,k)*(u_up(i,k)+weight(view)*CalMulti(V1,Vstar,k))/(u_down(i,k)+weight(view)*sum(U1(:,k))*CalMulti(V1,V1,k));
                    end
                    %Normalize U1 and V1 
                    v(k)=sum(U1(:,k));
                end
                Q1=diag(v);
                U1=U1/Q1;
                V1=V1*Q1;
                %Fixing Vstar ad U1, update V1
                v_up=X1'*U1;
                v_down=V1*(U1')*U1;
                for k=1:K
                    for j=1:viewNum
                        V1(j,k)=V1(j,k)*(v_up(j,k)+weight(view)*Vstar(j,k))/(v_down(j,k)+weight(view)*V1(j,k));
                    end
                end
                %O1_old=O1;
                O1=CalDif(X1,U1,V1,Q1,Vstar,weight(1));
                %plot_o1(it1)=O1;
                %plot_oo1(it1)=O1_old-O1;
                %fprintf('O1_old-O1=%f\n',O1_old-O1)
            end
            %O1
            %画图
        end
        %X2
        if view==2
            %while (abs(O2_old-O2)>1e-3)
            for it2=1:iter
                %Fixing Vstar&V2, update U2
                u_up=X2*V2;
                u_down=U2*(V2')*V2;
                for k=1:K
                    for i=1:cellNum2
                        U2(i,k)=U2(i,k)*(u_up(i,k)+weight(view)*CalMulti(V2,Vstar,k))/(u_down(i,k)+weight(view)*sum(U2(:,k))*CalMulti(V2,V2,k));
                    end
                %Normalize U2 and V2 
                v(k)=sum(U2(:,k));
                end
                Q2=diag(v);
                U2=U2/Q2;
                V2=V2*Q2;
                %Fixing Vstar ad U1, update V1
                v_up=(X2')*U2;
                v_down=V2*(U2')*U2;
                for k=1:K
                    for j=1:viewNum
                        V2(j,k)=V2(j,k)*(v_up(j,k)+weight(view)*Vstar(j,k))/(v_down(j,k)+weight(view)*V2(j,k));
                    end
                end
                %O2_old=O2;
                O2=CalDif(X2,U2,V2,Q2,Vstar,weight(2));
                %fprintf('O2_old-O2=%f\n',O2_old-O2)
            end
            %O2
        end
    end
        Vstar=weight(1)*V1*Q1+weight(2)*V2*Q2;
        O_old=O;
        O=CalDif(X2,U2,V2,Q2,Vstar,0.5)+CalDif(X1,U1,V1,Q1,Vstar,0.5);
        plot_o(iter)=O;
        fprintf('O_old-O=%f\n',O_old-O)
        fprintf('O=%f\n',O)
        
end
x=1:1:30;
plot(x,plot_o,'-*b')
axis([0,30,0,20]);
set(gca,'XTick',[0:5:30]);
set(gca,'YTick',[0:1:20]);
xlabel('迭代次数');
ylabel('误差值');
    
