%%Computing MECC multi-objective optimal designs in Application 2 (Table 3)
%%Dose-response models 

clear;
runningtime=cputime;  %record computation time
    tol = 10^(-3);
    N=201;         %number of design points   
    a=  -6.91;  b=  6.91;   %[a, b] is the design space
    e2=0.73;  e3=0.76;       %Efficiency used in the constraints
    
    u=linspace(a,b,N); %equally spaced N points in [a,b]                             
    v1=1.563; v2=0.825; v3=0.653; v4=0.137;  %true parameter values of theta 
    
    %Vectors and matrices used in the information matrices below. 
    p1=4;  %# of parameters in logistic model
    F1i=zeros(p1,p1,N);  
    for j=1:N 
        d=u(:,j);
        g=exp(v2*d+v3);
        f1=[1/(1+g) -v1*d*g/((1+g)^2)  -v1*g/((1+g)^2)  1];
        F1i(:,:,j)=(f1'*f1);
    end
   
%Compute the D-optimal design 
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A1(p1,p1); 
            for j=1:N
                A1 = A1+F1i(:,:,j)*w(j);
            end                   
            minimize( - det_rootn(A1) )   %D-optimality   
            0 <= w <= 1;
            sum(w)==1;
        cvx_end
        
        Loss1=1/(det(A1))^(1/p1);
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))'] %D-optimal design  
    
%Compute the ED50-optimal design 
        c1=[0, v3/(v2^2), -1/v2, 0]';
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A1(p1,p1); 
            for j=1:N
                A1 = A1+F1i(:,:,j)*w(j);
            end        
            minimize( matrix_frac(c1,A1) )   %ED50-optimality   
            0 <= w <= 1;
            sum(w)==1; 
        cvx_end

        Loss2=c1'*inv(A1)*c1;
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))'] %ED50-optimal design  
    
    %Compute the MED-optimal design 
        c2=[-1/(v1-1)/v2, (v3+log(v1-1))/v2^2, -1/v2, 0]';
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A1(p1,p1); 
            for j=1:N
                A1 = A1+F1i(:,:,j)*w(j);
            end
            minimize( matrix_frac(c2,A1) )   %MED-optimality   
            0 <= w <= 1;
            sum(w)==1; 
        cvx_end

        Loss3=c2'*inv(A1)*c2;
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))'] %MED-optimal design  
       
    %Compute the MECC efficiency constrained optimal design       
    cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A1(p1,p1); 
            for j=1:N
                A1 = A1+F1i(:,:,j)*w(j);
            end
            minimize( -det_rootn(A1) )   %D-optimality   
            0 <= w <= 1;
            sum(w)==1; 
            matrix_frac(c1,A1)-Loss2/e2 <=0;  %efficiency constraint
            matrix_frac(c2,A1)-Loss3/e3 <=0;  %efficiency constraint
        cvx_end

        Loss1d=1/((det(A1))^(1/p1));
        Loss2d=c1'*inv(A1)*c1;
        Loss3d=c2'*inv(A1)*c2;
        eff1=Loss1/Loss1d    %Efficiency at the MECC optimal design
        eff2=Loss2/Loss2d
        eff3=Loss3/Loss3d
        
    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))'] %optimal design  
    resulttime=cputime-runningtime  %computation time
 