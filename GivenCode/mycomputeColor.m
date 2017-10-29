function img = mycomputeColor(u,v)

% img = mycomputeColor(u,v) 
%
% input:
% u - first component of the flow field
% v - second component of the flow field
%
% output:
% img - rgb image representing the flow field

% saturation and value of the depiction are given by the size of the flow
% field; sizes are scaled to values between 0 and 1.
saturation = sqrt(u.^2+v.^2);
saturation_max = max(max(saturation));
saturation_scaled = saturation/saturation_max;

% hue is given by the direction of the flow field. The components of the
% flow field are interpreted as complex numbers (u + iv). As a first step,
% we compute their (principal) square root us + i vs.
us = sqrt((u+sqrt(u.^2+v.^2))/2);
vs = sign(v).*sqrt((-u+sqrt(u.^2+v.^2))/2);
% Now we define the hue as the argument of us + i vs, scaled to values
% between 0 and 1
hue = real((pi/2-atan(us./vs))/pi);
hue(hue==inf) = 0;
hue(hue==-inf) = 1;
hue(isnan(hue)) = 0.5;

% set up the flow field as hsv image
img_hsv = zeros([size(u),3]);
img_hsv(:,:,1) = hue;
img_hsv(:,:,2) = saturation_scaled;
img_hsv(:,:,3) = 1-(saturation_scaled.*(1-saturation_scaled)).^2;

% convert the hsv image into rgb
img = hsv2rgb(img_hsv);