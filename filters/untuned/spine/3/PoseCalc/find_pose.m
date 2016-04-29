% function find_pose
global object_points;
global A1;
global A2;
global A3;
global A4;
global Tcw1;
global Tcw2;
global Tcw3;
global Tcw4;
global Rcw1;
global Rcw2;
global Rcw3;
global Rcw4;
global measured;

red_blob = [-144;255;0];    % x;y;z position in mm w.r.t. body
green_blob = [-130;-265;0];
blue_blob = [315;0;0];
object_points = [red_blob, green_blob, blue_blob];



[fx1,fy1,cx1,cy1,Rcw1,Tcw1] = extract_camera_param(fd_client1);
[fx2,fy2,cx2,cy2,Rcw2,Tcw2] = extract_camera_param(fd_client2);
[fx3,fy3,cx3,cy3,Rcw3,Tcw3] = extract_camera_param(fd_client3);
[fx4,fy4,cx4,cy4,Rcw4,Tcw4] = extract_camera_param(fd_client4);
A_intrinsic_1 = [fx1 0 cx1; 0 fy1 cy1;0 0 1];
A_intrinsic_2 = [fx2 0 cx2; 0 fy2 cy2;0 0 1];
A_intrinsic_3 = [fx3 0 cx3; 0 fy3 cy3;0 0 1];
A_intrinsic_4 = [fx4 0 cx4; 0 fy4 cy4;0 0 1];
cam_param_1 = camera_parameters_store(A_intrinsic_1,Rcw1,Tcw1);
cam_param_2 = camera_parameters_store(A_intrinsic_2,Rcw2,Tcw2);
cam_param_3 = camera_parameters_store(A_intrinsic_3,Rcw3,Tcw3);
cam_param_4 = camera_parameters_store(A_intrinsic_4,Rcw4,Tcw4);

% load  first line of data
[tl1,tl2,tl3,tl4] = getlines(fd_client1,fd_client2,fd_client3,fd_client4);
s1 = extract_image_data(tl1,cam_param_1);
s2 = extract_image_data(tl2,cam_param_2);
s3 = extract_image_data(tl3,cam_param_3);
s4 = extract_image_data(tl4,cam_param_4);
z = [s1,s2,s3,s4];

A1 = A_intrinsic_1*1000; % convert mapping from pixel/mm to pixel/m
A2 = A_intrinsic_2*1000; 
A3 = A_intrinsic_3*1000; 
A4 = A_intrinsic_4*1000; 
Tcw1 = Tcw1/1000; % convert mapping from pixel/mm to pixel/m
Tcw2 = Tcw2/1000;
Tcw3 = Tcw3/1000;
Tcw4 = Tcw4/1000;

cam1_red_meas = [s1.red_data(1,1);s1.red_data(1,2)];
cam1_green_meas = [s1.green_data(1,1);s1.green_data(1,2)];
cam1_blue_meas = [s1.blue_data(1,1);s1.blue_data(1,2)];
cam2_red_meas = [s2.red_data(1,1);s2.red_data(1,2)];
cam2_green_meas = [s2.green_data(1,1);s2.green_data(1,2)];
cam2_blue_meas = [s2.blue_data(1,1);s2.blue_data(1,2)];
cam3_red_meas = [s3.red_data(1,1);s3.red_data(1,2)];
cam3_green_meas = [s3.green_data(1,1);s3.green_data(1,2)];
cam3_blue_meas = [s3.blue_data(1,1);s3.blue_data(1,2)];
cam4_red_meas = [s4.red_data(1,1);s4.red_data(1,2)];
cam4_green_meas = [s4.green_data(1,1);s4.green_data(1,2)];
cam4_blue_meas = [s4.blue_data(1,1);s4.blue_data(1,2)];
measured = [cam1_red_meas;cam1_green_meas;cam1_blue_meas;
    cam2_red_meas;cam2_green_meas;cam2_blue_meas;
    cam3_red_meas;cam3_green_meas;cam3_blue_meas;
    cam4_red_meas;cam4_green_meas;cam4_blue_meas];


x0 = [0;0;0;0;0;0;0;0;0;]; % initial guess

[x,resnorm] = lsqnonlin(@pixelError,x0);
for i = 1:sampleLength
% while ~feof(fd_client1)
    
    % load next line
    [tl1,tl2,tl3,tl4] = getlines(fd_client1,fd_client2,fd_client3,fd_client4);
    s1 = extract_image_data(tl1,cam_param_1);
    s2 = extract_image_data(tl2,cam_param_2);
    s3 = extract_image_data(tl3,cam_param_3);
    s4 = extract_image_data(tl4,cam_param_4);
    
    if (s1.number_red == 0)
        cam1_red_meas = [-1;-1];
    else
        cam1_red_meas = [s1.red_data(1,1);s1.red_data(1,2)];
    end
    if (s1.number_green == 0)
        cam1_green_meas = [-1;-1];
    else
        cam1_green_meas = [s1.green_data(1,1);s1.green_data(1,2)];
    end
    if (s1.number_blue == 0)
        cam1_blue_meas = [-1;-1];
    else
        cam1_blue_meas = [s1.blue_data(1,1);s1.blue_data(1,2)];
    end
    
    if (s2.number_red == 0)
        cam2_red_meas = [-1;-1];
    else
        cam2_red_meas = [s2.red_data(1,1);s2.red_data(1,2)];
    end
    if (s2.number_green == 0)
        cam2_green_meas = [-1;-1];
    else
        cam2_green_meas = [s2.green_data(1,1);s2.green_data(1,2)];
    end
    if (s2.number_blue == 0)
        cam2_blue_meas = [-1;-1];
    else
        cam2_blue_meas = [s2.blue_data(1,1);s2.blue_data(1,2)];
    end    
    
    if (s3.number_red == 0)
        cam3_red_meas = [-1;-1];
    else
        cam3_red_meas = [s3.red_data(1,1);s3.red_data(1,2)];
    end
    if (s3.number_green == 0)
        cam3_green_meas = [-1;-1];
    else
        cam3_green_meas = [s3.green_data(1,1);s3.green_data(1,2)];
    end
    if (s3.number_blue == 0)
        cam3_blue_meas = [-1;-1];
    else
        cam3_blue_meas = [s3.blue_data(1,1);s3.blue_data(1,2)];
    end    
        
    if (s4.number_red == 0)
        cam4_red_meas = [-1;-1];
    else
        cam4_red_meas = [s4.red_data(1,1);s4.red_data(1,2)];
    end
    if (s4.number_green == 0)
        cam4_green_meas = [-1;-1];
    else
        cam4_green_meas = [s4.green_data(1,1);s4.green_data(1,2)];
    end
    if (s4.number_blue == 0)
        cam4_blue_meas = [-1;-1];
    else
        cam4_blue_meas = [s4.blue_data(1,1);s4.blue_data(1,2)];
    end    
    
    measured = [cam1_red_meas;cam1_green_meas;cam1_blue_meas;
        cam2_red_meas;cam2_green_meas;cam2_blue_meas;
        cam3_red_meas;cam3_green_meas;cam3_blue_meas;
        cam4_red_meas;cam4_green_meas;cam4_blue_meas];
    
    [x,resnorm] = lsqnonlin(@pixelError,x);
    out(:,i) = x;
end









