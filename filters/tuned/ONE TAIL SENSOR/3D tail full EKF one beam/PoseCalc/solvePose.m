function x = solvePose(x0,z,object_points)
    x(:,1) = x0;
    for i = 1:15
        [x_est,norm_err] = solvePoseEst(x(:,i),z,object_points)   
        x(:,i+1) = x_est;
    end
end