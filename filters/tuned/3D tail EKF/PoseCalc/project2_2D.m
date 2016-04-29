function [x] = project2_2D(A,Rcw,Tcw,Twp)
    image_homog = A*Rcw*(Twp - Tcw);
    x(1,1) = image_homog(1)/image_homog(3);
    x(2,1) = image_homog(2)/image_homog(3);
end