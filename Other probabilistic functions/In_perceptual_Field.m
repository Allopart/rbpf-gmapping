function output= In_perceptual_Field(ii,jj,pose,y)

xc = ii-0.5;
yc = jj-0.5; %There will probably be a problem with this line
%x=pose(1), y=pose(2), theta=pose(3)
r=sqrt((xc-pose(1))^2+(yc-pose(2))^2);
theta=atan2(yc-pose(2), xc-pose(1))-pose(3);


if ((theta > max(y(1,:))) || (theta < min(y(1,:))) || (r > max(y(2,:))))
    output=0;
else
    output=1;
end

