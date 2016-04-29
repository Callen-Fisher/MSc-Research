function y = pixelError(x)
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
    
    R = x(1:3);
    G = x(4:6);
    B = x(7:9);
    cam1_red = project2_2D(A1,Rcw1,Tcw1,R);
    cam1_green = project2_2D(A1,Rcw1,Tcw1,G);
    cam1_blue = project2_2D(A1,Rcw1,Tcw1,B);
    cam2_red = project2_2D(A2,Rcw2,Tcw2,R);
    cam2_green = project2_2D(A2,Rcw2,Tcw2,G);
    cam2_blue = project2_2D(A2,Rcw2,Tcw2,B);
    cam3_red = project2_2D(A3,Rcw3,Tcw3,R);
    cam3_green = project2_2D(A3,Rcw3,Tcw3,G);
    cam3_blue = project2_2D(A3,Rcw3,Tcw3,B);
    cam4_red = project2_2D(A4,Rcw4,Tcw4,R);
    cam4_green = project2_2D(A4,Rcw4,Tcw4,G);
    cam4_blue = project2_2D(A4,Rcw4,Tcw4,B);    
    
    estimated = [cam1_red;cam1_green;cam1_blue;
    cam2_red;cam2_green;cam2_blue;
    cam3_red;cam3_green;cam3_blue;
    cam4_red;cam4_green;cam4_blue];
    
    n = length(measured);
    k = 1;    
    for i = 1:n
        if (measured(i,1) ~= -1)
            measured_mod(k,1) = measured(i,1);
            estimated_mod(k,1) = estimated(i,1);
            k = k + 1;
        end
    end
    y = measured_mod - estimated_mod;    
end