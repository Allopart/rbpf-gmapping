function result = gpc_test(real_x, data)
% Tests the gpc code

	
	p = [1 0; 0 1; -1 0; 2 1; 4 2]';
	alpha = [deg2rad(0);deg2rad(10);deg2rad(20);deg2rad(50);deg2rad(-20);];

	noise = 0;
	for i=1:size(data,2)
		corr{i}.p = data(:,i);
		corr{i}.q = rot(theta)*p(:,i)+t + randn(2,1)*noise;
		corr{i}.C = vers(alpha(i))*vers(alpha(i))';
	%	corr{i}.C = eye(2);
	end

	res = gpc(corr);
    result=res.x;
% 	fprintf('Real     : %s\n', sprintf('%f ', real_x));
% 	fprintf('Estimated: %s\n', sprintf('%f ', res.x));


function R = rot(theta)
	R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
	
function v = vers(theta)
	v = [cos(theta); sin(theta)];