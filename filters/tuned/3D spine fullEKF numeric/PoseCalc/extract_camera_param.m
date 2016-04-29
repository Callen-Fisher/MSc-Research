function [fx,fy,cx,cy,Rcw,Tcw] = extract_camera_param(fd)
    % get intrinsics
    tln = fgetl(fd);
    [token,remain] = strtok(tln, ':');
    [token,remain] = strtok(remain, ' ');
    [token,remain] = strtok(remain, ',');
    fx = str2num(token);
    [token,remain] = strtok(remain, ',');
    fy = str2num(token);
    [token,remain] = strtok(remain, ',');
    cx = str2num(token);
    [token,remain] = strtok(remain, ',');
    cy = str2num(token);
    
    %get rotation matrix
    tln = fgetl(fd);
    [token,remain] = strtok(tln, ':');
    [token,remain] = strtok(remain, ' ');
    [token,remain] = strtok(remain, ',');
    theta_z = str2num(token);
    [token,remain] = strtok(remain, ',');
    theta_x = str2num(token);
    [token,remain] = strtok(remain, ',');
    theta_z2 = str2num(token);
    
    Rz = [cos(theta_z) sin(theta_z) 0;
          -sin(theta_z) cos(theta_z) 0;
          0 0 1];
    Rx = [1 0 0;
          0 cos(theta_x) sin(theta_x);
          0 -sin(theta_x) cos(theta_x);];
    Rz2 = [cos(theta_z2) sin(theta_z2) 0;
          -sin(theta_z2) cos(theta_z2) 0;
          0 0 1];      

    Rcw = Rz2*Rx*Rz;
    
    %get translation matrix
    tln = fgetl(fd);
    [token,remain] = strtok(tln, ':');
    [token,remain] = strtok(remain, ' ');
    [token,remain] = strtok(remain, ',');
    Tcw_x = str2num(token);
    [token,remain] = strtok(remain, ',');
    Tcw_y = str2num(token);
    [token,remain] = strtok(remain, ',');
    Tcw_z = str2num(token);
    Tcw = [Tcw_x;Tcw_y;Tcw_z];
end