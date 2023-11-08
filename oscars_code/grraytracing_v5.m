% This is mostly my own. 
% For the one inspired purely by the paper, see v1
% For the one thats paper + my own geodesics, see v2
    % Couldn't get that one to work. Let's see if I can get this one
% For CPU acceleration, see v4 or v5

clear

% Camera specifications- leaving out time component for now...
camPos = [-25 0 3]*0.7;
targetPos = [0 0 0];

FOV = pi/2;
verticalPixels = 100;
horizontalPixels = 100;

%Formulate vectors used for camera transformation
vertical = [0 -2.54 6];
camDir = targetPos - camPos;
horizontal = cross(vertical,camDir);

unitCamDir = camDir./norm(camDir);
unitHorizontal = horizontal./norm(horizontal);
unitVertical = cross(unitCamDir,unitHorizontal);

outerPixelHorzDist = tan(FOV/2);
outerPixelVertDist = outerPixelHorzDist*(verticalPixels-1)/(horizontalPixels-1);

horizontalStep = 2*outerPixelHorzDist/(horizontalPixels-1) .* unitHorizontal;
verticalStep = 2*outerPixelVertDist/(verticalPixels-1) .* unitVertical;

ray_11 = unitCamDir - outerPixelHorzDist.*unitHorizontal + outerPixelVertDist.*unitVertical;

%Define system of diffeqs
metricName = 'Schwarz';
% metricName = 'Kerr';

% load(['geodesics_' metricName '.mat']);
load(['geodesics_lagrange_' metricName '.mat'])

eq1 = geodesics{1};
eq2 = geodesics{2};
eq3 = geodesics{3};
eq4 = geodesics{4};

eqs = [eq1; eq2; eq3; eq4];

[V,S] = odeToVectorField(eqs);

geoFunc = matlabFunction(V,'vars',{'lambda','Y'});

%Drawing the screen
viewPort = zeros(verticalPixels*horizontalPixels,3);
diskPic = imread('diskTexture_3.png');
background = imread('stars.jpg');
% background = imread('stars2.tif');
dimsBG = size(background);
dimsDisk = size(diskPic);

lambda_f = 100;

diskRadius = r_s + 15;

tic
parfor pixel = 1:horizontalPixels*verticalPixels
    warning('off','MATLAB:ode23t:IntegrationTolNotMet');
    i = mod(pixel-1,horizontalPixels)+1;
    j = ceil(pixel/verticalPixels);

    %Ray direction
    thisRay = ray_11 + horizontalStep*j - verticalStep*i;
    
    thisRaySph = cartToSphereVel(thisRay,camPos);
    camPosSph = cartToSpherePos(camPos);

    timeComp = norm(thisRay);

    %Define initial conditions based on camera pos and ray dir
    initCond = [camPosSph(1); thisRaySph(1); 0; timeComp; ...
        camPosSph(2); thisRaySph(2); camPosSph(3); thisRaySph(3)];
    
    soln = ode23t(geoFunc,[0 lambda_f], initCond);

    theseRadii = soln.y(1,:);
    radMask = (theseRadii < diskRadius) & (theseRadii > (r_s+1));
    
    theseThetas = soln.y(5,:);
    delThetas = [0 theseThetas-pi/2] .* [theseThetas-pi/2 0];
    thetaMask = delThetas < 0;

    thesePhis = soln.y(7,:);

    accretionMask = radMask & thetaMask(2:end);
%     accretionMask = 0;

    if soln.x(end) < (lambda_f - 10)
        backColor = zeros(1,1,3);
    else
        thisTheta = soln.y(5,end);
        thisPhi = soln.y(7,end);

        beta = 0.400458119246;
        a = [sin(thisTheta)*cos(thisPhi), sin(thisTheta)*sin(thisPhi), cos(thisTheta)];
        b = cos(beta)*a + sin(beta)*[sin(thisTheta)*cos(thisPhi)*(1-cos(beta)),a(3),-a(2)];

        u = atan(sqrt(b(1).^2 + b(2).^2) ./ b(3));
        v = atan(b(2)/b(1));

        raw_u = u/pi;
        raw_v = 1 - (v/(2*pi)+0.5);

        u = mod(floor(raw_u*dimsBG(1)),dimsBG(1))+1;
        v = mod(floor(raw_v*dimsBG(2)),dimsBG(2))+1;

        backColor = double(background(u,v,:)) * 0.4;
    end

    if any(accretionMask)
        rIntersections = theseRadii(accretionMask);
        thisR = rIntersections(1);
        phiIntersections = thesePhis(accretionMask);
        thisPhi = phiIntersections(1);
        
        normalizedR = (thisR/diskRadius);

        raw_u = normalizedR/2*cos(thisPhi)+0.5;
        raw_v = normalizedR/2*sin(thisPhi)+0.5;

        u = mod(floor(raw_u*dimsDisk(1)),dimsDisk(1))+1;
        v = mod(floor(raw_v*dimsDisk(2)),dimsDisk(2))+1;

        accretionColor = double(diskPic(u,v,:));
        
        opacity = exp(1-1/(1-normalizedR^2));

        if sum(accretionMask) > 1
            thisR_2 = rIntersections(2);
            thisPhi_2 = phiIntersections(2);

            normalizedR_2 = (thisR_2/diskRadius);

            raw_u_2 = normalizedR_2/2*cos(thisPhi_2)+0.5;
            raw_v_2 = normalizedR_2/2*sin(thisPhi_2)+0.5;

            u_2 = mod(floor(raw_u_2*dimsDisk(1)),dimsDisk(1))+1;
            v_2 = mod(floor(raw_v_2*dimsDisk(2)),dimsDisk(2))+1;

            accretionColor_2 = double(diskPic(u_2,v_2,:)) * 2;

            opacity_2 = exp(1-1/(1-normalizedR_2^2));

            thisColor_2 = accretionColor_2*opacity_2 + backColor*(1-opacity_2);
            thisColor = accretionColor*opacity + thisColor_2*(1-opacity);
        else
            thisColor = accretionColor*opacity + backColor*(1-opacity);
        end
    else
        thisColor = backColor;
    end

    viewPort(pixel,:) = thisColor;
end
toc

viewPort = reshape(viewPort,horizontalPixels,verticalPixels,3);
image(uint8(viewPort));

function [spherical] = cartToSphereVel(cart,pos)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    r = sqrt(x^2+y^2+z^2);
    J = [x/r y/r z/r; ...
        x*z/(r^2*sqrt(x^2+y^2)) y*z/(r^2*sqrt(x^2+y^2)) -(x^2+y^2)/(r^2*sqrt(x^2+y^2)); ...
        -y/(x^2+y^2) x/(x^2+y^2) 0];
    spherical = J*cart.';
    spherical = spherical.';
end

function [spherical] = cartToSpherePos(cart)
    x = cart(1);
    y = cart(2);
    z = cart(3);
    r = sqrt(x^2+y^2+z^2);
    theta = acos(z/r);
    if x > 0
        phi = atan(y/x);
    elseif x < 0 && y >= 0
        phi = atan(y/x) + pi;
    elseif x < 0 && y < 0
        phi = atan(y/x) - pi;
    elseif x == 0 && y > 0
        phi = pi/2;
    elseif x == 0 && y < 0
        phi = -pi/2;
    elseif x == 0 && y == 0
        phi = 0;
    end
    spherical = [r theta phi];
end