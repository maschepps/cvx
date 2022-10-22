%Compute maximin D-optimal design in Application 3
%%Four dose-response models 

clear;
runningtime=cputime;  %record computation time

    tol = 10^(-4);
    N=130;         %number of design points   
    a=0;  b=  500; %[a, b] is the design space
    u=linspace(a,b,N); %equally spaced N points in [a,b]
    v0=60; v1=294; v2=25;   %true parameter values for Emax I model 
    v02=60;  v12=340; v22=107.14; %true parameter values for Emax II model 
    v03=49.62; v13=290.51; v23=150; v33=45.51; %parameter values for logistic model
  
    %Vectors and matrices are used in the information matrices below.
    p1=2; p2=3; p3=p2; p4=4; %# of parameters in the four models 
    F1i=zeros(p1,p1,N); F2i=zeros(p2,p2,N); F3i=zeros(p3,p3,N); F4i=zeros(p4,p4,N);  
    for j=1:N 
        d=u(:,j);
        f1=[1 d];
        f2=[1 d/(v2+d) -v1*d/(v2+d)^2];
        f3=[1 d/(v22+d) -v12*d/(v22+d)^2];
        g=exp((v23-d)/v33);
        f4=[1 1/(1+g) -v13*g/v33/((1+g)^2)  v13*g*(v23-d)/v33^2/((1+g)^2)];
        F1i(:,:,j)=(f1'*f1);
        F2i(:,:,j)=(f2'*f2);
        F3i(:,:,j)=(f3'*f3);
        F4i(:,:,j)=(f4'*f4);
    end
  
%Compute the D-optimal design for linear model
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A1(p1,p1); 
            for j=1:N
                A1 = A1+F1i(:,:,j)*w(j);
            end        
            minimize( - det_rootn(A1) )      
            0 <= w <= 1;
            sum(w)==1;
        cvx_end
        Loss1=1/(det(A1))^(1/p1);
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   
    
%Compute the D-optimal design for Emax I model
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A2(p2,p2); 
            for j=1:N
                A2 = A2+F2i(:,:,j)*w(j);
            end        
            minimize( - det_rootn(A2) )   %D-optimality   
            0 <= w <= 1;
            sum(w)==1;
        cvx_end
        Loss2=1/(det(A2))^(1/p2);
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   
     
    %Compute the D-optimal design for Emax II model
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A3(p3,p3); 

            for j=1:N
                A3 = A3+F3i(:,:,j)*w(j);
            end               
            minimize( - det_rootn(A3) )   %D-optimality   
            0 <= w <= 1;
            sum(w)==1;
        cvx_end
        Loss3=1/(det(A3))^(1/p3);
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   
    
    %Compute the D-optimal design for logistic model
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A4(p4,p4); 
            for j=1:N
                A4 = A4+F4i(:,:,j)*w(j);
            end      
            minimize( - det_rootn(A4) )   %D-optimality   
            0 <= w <= 1;
            sum(w)==1;
        cvx_end
        Loss4=1/(det(A4))^(1/p4);
     %design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   
    
     %Compute the Maximin D-efficiency design for the 4 models
        cvx_begin
            cvx_precision high
            variable w(1,N+1);
            expression B1(p1,p1); 
            expression B2(p2,p2);
            expression B3(p3,p3); 
            expression B4(p4,p4); 
            for j=1:N
                B1 = B1+F1i(:,:,j)*w(j);
                B2 = B2+F2i(:,:,j)*w(j);
                B3 = B3+F3i(:,:,j)*w(j);
                B4 = B4+F4i(:,:,j)*w(j);
            end         
            maximize w(N+1)
            det_rootn(B1)- w(N+1)/Loss1 >=0; 
            det_rootn(B2)- w(N+1)/Loss2 >=0; 
            det_rootn(B3)- w(N+1)/Loss3 >=0; 
            det_rootn(B4)- w(N+1)/Loss4 >=0; 
            0 <= w;
            sum(w(1:N))==1;
        cvx_end

        Loss1d=1/(det(B1))^(1/p1);
        Loss2d=1/(det(B2))^(1/p2);
        Loss3d=1/(det(B3))^(1/p3);
        Loss4d=1/(det(B4))^(1/p4);
        eff1=Loss1/Loss1d  %Efficiency at the maximin design
        eff2=Loss2/Loss2d
        eff3=Loss3/Loss3d
        eff4=Loss4/Loss4d
        w=w(1:N);
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))'] %Maximin design   
    
   resulttime=cputime-runningtime  %computation time
 
    