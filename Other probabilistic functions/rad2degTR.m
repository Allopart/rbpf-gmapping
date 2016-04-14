function rotTR= rad2degTR(TR)
    rotTR=zeros(size(TR));
    rotTR(1,1)=rad2deg(acos(TR(1,1)));
    rotTR(1,2)=-rad2deg(asin(TR(1,2)));
    rotTR(2,1)=rad2deg(asin(TR(2,1)));
    rotTR(2,2)=rad2deg(acos(TR(2,2)));
end