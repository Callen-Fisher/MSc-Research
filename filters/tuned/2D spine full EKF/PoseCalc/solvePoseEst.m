function [x_est,resnorm] = solvePoseEst(x,z,object_points)    
    num_measurements = 0;
    % only allow single blob detection for now
    for i = 1:4
        if (z(i).number_red == 1)
            num_measurements = num_measurements + 2;
        end
        if (z(i).number_green == 1)
            num_measurements = num_measurements + 2;
        end
        if (z(i).number_blue == 1)
            num_measurements = num_measurements + 2;
        end        
    end
    
    A_lin = zeros(num_measurements,7);    
    b_lin = zeros(num_measurements,1);
    
    z_est = zeros(num_measurements,1);
    z_meas = zeros(num_measurements,1);    
    
    Twp = x(1:3);
    qwp = x(4:7);
    k = 1;
    for i = 1:4
        A = z(1,i).cam_param.intrinsic_mat*1000; % convert mapping from pixel/mm to pixel/m
        Rcw = z(1,i).cam_param.rotation_w2c;
        Tcw = z(1,i).cam_param.translation_w2c/1000; % convert mapping from pixel/mm to pixel/m
        % red blobs
        if (z(1,i).number_red == 1)
            % do measurement estimation
            u = project2_2D(A,Rcw,Tcw,Twp,qwp,object_points(:,1)/1000);
            z_est(k:k+1,1) = u;            
            
            % do jacobian calculation
            K = A*Rcw;
            X = object_points(1,1)/1000;
            Y = object_points(2,1)/1000;
            Z = object_points(3,1)/1000;
            
            % calculate B's
            B1 = (qwp(1)^2+qwp(2)^2-qwp(3)^2-qwp(4)^2)*X + ...
                + 2*(qwp(2)*qwp(3) - qwp(1)*qwp(4))*Y + ...
                + 2*(qwp(2)*qwp(4) + qwp(1)*qwp(3))*Z + ...
                + Twp(1) - z(1,i).cam_param.translation_w2c(1)/1000; % convert mapping from pixel/mm to pixel/m
            B2 = 2*(qwp(2)*qwp(3) + qwp(1)*qwp(4))*X + ...
                + (qwp(1)^2-qwp(2)^2+qwp(3)^2-qwp(4)^2)*Y + ...
                + 2*(qwp(3)*qwp(4) - qwp(1)*qwp(2))*Z + ...
                + Twp(2) - z(1,i).cam_param.translation_w2c(2)/1000; % convert mapping from pixel/mm to pixel/m
            B3 = 2*(qwp(2)*qwp(4) - qwp(1)*qwp(3))*X + ...
                + 2*(qwp(3)*qwp(4) + qwp(1)*qwp(2))*Y + ...
                + (qwp(1)^2-qwp(2)^2-qwp(3)^2+qwp(4)^2)*Z + ...
                + Twp(3) - z(1,i).cam_param.translation_w2c(3)/1000; % convert mapping from pixel/mm to pixel/m            
            
            b = K(3,1)*B1 + K(3,2)*B2 + K(3,3)*B3;
            bu = b;
            bv = b;
            au = K(1,1)*B1 + K(1,2)*B2 + K(1,3)*B3;
            av = K(2,1)*B1 + K(2,2)*B2 + K(2,3)*B3;
            
            dB1_q = [2*(qwp(1)*X - qwp(4)*Y + qwp(3)*Z),2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z),...
                    2*(-qwp(3)*X + qwp(2)*Y + qwp(1)*Z),2*(-qwp(4)*X - qwp(1)*Y + qwp(2)*Z)];
                
            dB2_q = [2*(qwp(4)*X + qwp(1)*Y - qwp(2)*Z),2*(qwp(3)*X - qwp(2)*Y - qwp(1)*Z),...
                    2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z),2*(qwp(1)*X - qwp(4)*Y + qwp(3)*Z)];
                
            dB3_q = [2*(-qwp(3)*X + qwp(2)*Y + qwp(1)*Z),2*(qwp(4)*X + qwp(1)*Y - qwp(2)*Z),...
                    2*(-qwp(1)*X + qwp(4)*Y - qwp(3)*Z),2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z)];
                
            dB1_x = [1 0 0 dB1_q];
            dB2_x = [0 1 0 dB2_q];
            dB3_x = [0 0 1 dB3_q];
            
            au_prime = K(1,1)*dB1_x + K(1,2)*dB2_x + K(1,3)*dB3_x;
            av_prime = K(2,1)*dB1_x + K(2,2)*dB2_x + K(2,3)*dB3_x;
            b_prime = K(3,1)*dB1_x + K(3,2)*dB2_x + K(3,3)*dB3_x;
            bu_prime = b_prime;
            bv_prime = b_prime;            
            
            b_sq = b^2;
            du_dx = (au_prime*bu - au*bu_prime)/b_sq;
            dv_dx = (av_prime*bv - av*bv_prime)/b_sq;
            A_lin(k:k+1,:) = [du_dx;dv_dx];     
            
            z_meas(k,1) = z(1,i).red_data(1,1);
            z_meas(k+1,1) = z(1,i).red_data(1,2);             
            
            k = k + 2;
        end
        % green blobs
        if (z(1,i).number_green == 1)
            % do measurement estimation
            u = project2_2D(A,Rcw,Tcw,Twp,qwp,object_points(:,2)/1000);
            z_est(k:k+1,1) = u;
            
            % do jacobian calc
            K = A*Rcw;
            X = object_points(1,2)/1000; % convert mapping from pixel/mm to pixel/m
            Y = object_points(2,2)/1000; % convert mapping from pixel/mm to pixel/m
            Z = object_points(3,2)/1000;
            
            % calculate B's
            B1 = (qwp(1)^2+qwp(2)^2-qwp(3)^2-qwp(4)^2)*X + ...
                + 2*(qwp(2)*qwp(3) - qwp(1)*qwp(4))*Y + ...
                + 2*(qwp(2)*qwp(4) + qwp(1)*qwp(3))*Z + ...
                + Twp(1) - z(1,i).cam_param.translation_w2c(1)/1000; % convert mapping from pixel/mm to pixel/m
            B2 = 2*(qwp(2)*qwp(3) + qwp(1)*qwp(4))*X + ...
                + (qwp(1)^2-qwp(2)^2+qwp(3)^2-qwp(4)^2)*Y + ...
                + 2*(qwp(3)*qwp(4) - qwp(1)*qwp(2))*Z + ...
                + Twp(2) - z(1,i).cam_param.translation_w2c(2)/1000; % convert mapping from pixel/mm to pixel/m
            B3 = 2*(qwp(2)*qwp(4) - qwp(1)*qwp(3))*X + ...
                + 2*(qwp(3)*qwp(4) + qwp(1)*qwp(2))*Y + ...
                + (qwp(1)^2-qwp(2)^2-qwp(3)^2+qwp(4)^2)*Z + ...
                + Twp(3) - z(1,i).cam_param.translation_w2c(3)/1000; % convert mapping from pixel/mm to pixel/m            
            
            b = K(3,1)*B1 + K(3,2)*B2 + K(3,3)*B3;
            bu = b;
            bv = b;
            au = K(1,1)*B1 + K(1,2)*B2 + K(1,3)*B3;
            av = K(2,1)*B1 + K(2,2)*B2 + K(2,3)*B3;
                
            dB1_q = [2*(qwp(1)*X - qwp(4)*Y + qwp(3)*Z),2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z),...
                    2*(-qwp(3)*X + qwp(2)*Y + qwp(1)*Z),2*(-qwp(4)*X - qwp(1)*Y + qwp(2)*Z)];
            dB2_q = [2*(qwp(4)*X + qwp(1)*Y - qwp(2)*Z),2*(qwp(3)*X - qwp(2)*Y - qwp(1)*Z),...
                    2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z),2*(qwp(1)*X - qwp(4)*Y + qwp(3)*Z)];
            dB3_q = [2*(-qwp(3)*X + qwp(2)*Y + qwp(1)*Z),2*(qwp(4)*X + qwp(1)*Y - qwp(2)*Z),...
                    2*(-qwp(1)*X + qwp(4)*Y - qwp(3)*Z),2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z)];
            
            dB1_x = [1 0 0 dB1_q];
            dB2_x = [0 1 0 dB2_q];
            dB3_x = [0 0 1 dB3_q];
            
            au_prime = K(1,1)*dB1_x + K(1,2)*dB2_x + K(1,3)*dB3_x;
            av_prime = K(2,1)*dB1_x + K(2,2)*dB2_x + K(2,3)*dB3_x;
            b_prime = K(3,1)*dB1_x + K(3,2)*dB2_x + K(3,3)*dB3_x;
            bu_prime = b_prime;
            bv_prime = b_prime;            
            
            b_sq = b^2;
            du_dx = (au_prime*bu - au*bu_prime)/b_sq;
            dv_dx = (av_prime*bv - av*bv_prime)/b_sq;
            A_lin(k:k+1,:) = [du_dx;dv_dx];
            
            z_meas(k,1) = z(1,i).green_data(1,1);
            z_meas(k+1,1) = z(1,i).green_data(1,2);            
            
            k = k + 2;
        end
        % blue blobs
        if (z(1,i).number_blue == 1)
            u = project2_2D(A,Rcw,Tcw,Twp,qwp,object_points(:,3)/1000);
            z_est(k:k+1,1) = u;
            
            K = A*Rcw;
            X = object_points(1,3)/1000; % convert mapping from pixel/mm to pixel/m
            Y = object_points(2,3)/1000; % convert mapping from pixel/mm to pixel/m
            Z = object_points(3,3)/1000;

            % calculate B's
            B1 = (qwp(1)^2+qwp(2)^2-qwp(3)^2-qwp(4)^2)*X + ...
                + 2*(qwp(2)*qwp(3) - qwp(1)*qwp(4))*Y + ...
                + 2*(qwp(2)*qwp(4) + qwp(1)*qwp(3))*Z + ...
                + Twp(1) - z(1,i).cam_param.translation_w2c(1)/1000; % convert mapping from pixel/mm to pixel/m
            B2 = 2*(qwp(2)*qwp(3) + qwp(1)*qwp(4))*X + ...
                + (qwp(1)^2-qwp(2)^2+qwp(3)^2-qwp(4)^2)*Y + ...
                + 2*(qwp(3)*qwp(4) - qwp(1)*qwp(2))*Z + ...
                + Twp(2) - z(1,i).cam_param.translation_w2c(2)/1000; % convert mapping from pixel/mm to pixel/m
            B3 = 2*(qwp(2)*qwp(4) - qwp(1)*qwp(3))*X + ...
                + 2*(qwp(3)*qwp(4) + qwp(1)*qwp(2))*Y + ...
                + (qwp(1)^2-qwp(2)^2-qwp(3)^2+qwp(4)^2)*Z + ...
                + Twp(3) - z(1,i).cam_param.translation_w2c(3)/1000; % convert mapping from pixel/mm to pixel/m            
            
            b = K(3,1)*B1 + K(3,2)*B2 + K(3,3)*B3;
            bu = b;
            bv = b;
            au = K(1,1)*B1 + K(1,2)*B2 + K(1,3)*B3;
            av = K(2,1)*B1 + K(2,2)*B2 + K(2,3)*B3;
                
            dB1_q = [2*(qwp(1)*X - qwp(4)*Y + qwp(3)*Z),2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z),...
                    2*(-qwp(3)*X + qwp(2)*Y + qwp(1)*Z),2*(-qwp(4)*X - qwp(1)*Y + qwp(2)*Z)];
            dB2_q = [2*(qwp(4)*X + qwp(1)*Y - qwp(2)*Z),2*(qwp(3)*X - qwp(2)*Y - qwp(1)*Z),...
                    2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z),2*(qwp(1)*X - qwp(4)*Y + qwp(3)*Z)];
            dB3_q = [2*(-qwp(3)*X + qwp(2)*Y + qwp(1)*Z),2*(qwp(4)*X + qwp(1)*Y - qwp(2)*Z),...
                    2*(-qwp(1)*X + qwp(4)*Y - qwp(3)*Z),2*(qwp(2)*X + qwp(3)*Y + qwp(4)*Z)];
            
            dB1_x = [1 0 0 dB1_q];
            dB2_x = [0 1 0 dB2_q];
            dB3_x = [0 0 1 dB3_q];
            
            au_prime = K(1,1)*dB1_x + K(1,2)*dB2_x + K(1,3)*dB3_x;
            av_prime = K(2,1)*dB1_x + K(2,2)*dB2_x + K(2,3)*dB3_x;
            b_prime = K(3,1)*dB1_x + K(3,2)*dB2_x + K(3,3)*dB3_x;
            bu_prime = b_prime;
            bv_prime = b_prime;
            
            b_sq = b^2;
            du_dx = (au_prime*bu - au*bu_prime)/b_sq;
            dv_dx = (av_prime*bv - av*bv_prime)/b_sq;
            A_lin(k:k+1,:) = [du_dx;dv_dx];
            
            z_meas(k,1) = z(1,i).blue_data(1,1);
            z_meas(k+1,1) = z(1,i).blue_data(1,2);                        

            k = k + 2;
        end
    end    
    
    b_lin = z_meas - z_est;
    
    [delta_x,norm_err,delta_x_plus_x] = solveLSQ(A_lin,b_lin,x);
    x_est = delta_x_plus_x;
    resnorm = norm_err;
end