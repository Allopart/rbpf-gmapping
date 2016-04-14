function x_new= Update_Pos (x_old, TR, TT)
%     x_new=zeros(3,1);
%     x_new(1)=TR(1,1)*x_old(1)+TR(1,2)*x_old(2)+TT(1);
%     x_new(2)=TR(2,1)*x_old(1)+TR(2,2)*x_old(2)+TT(2);
%     x_new(3)=x_old(3);
     
%     T=[TR TT; 0 0 1];
%     x_new=T*[x_old(1); x_old(2);1];
    
    x_new(1)=x_old(1)+TT(1);
    x_new(2)=x_old(2)+TT(2);
    theta=sin(TR(2,1));
    x_new(3)=x_old(3)+theta;
end