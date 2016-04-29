function s = camera_parameters_store(A,Rcw,Tcw)
    f1 = 'intrinsic_mat';
    f2 = 'rotation_w2c';
    f3 = 'translation_w2c';
    
    v1 = A;
    v2 = Rcw;
    v3 = Tcw;
    s = struct(f1,v1,f2,v2,f3,v3);
end