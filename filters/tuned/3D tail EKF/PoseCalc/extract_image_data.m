function s=extract_image_data(image_data,param)
    % get id
    [token,remain] = strtok(image_data, '*C:');
    id = str2num(token);
    
    % get number of red blobs
    [token,remain] = strtok(remain, ':*,|');
    r = str2num(token);
    % get number of green blobs
    [token,remain] = strtok(remain, ':*,|');
    g = str2num(token);
    % get number of blue blobs
    [token,remain] = strtok(remain, ':*,|');
    b = str2num(token);
    
    red_data = zeros(5,2);
    green_data = zeros(5,2);
    blue_data = zeros(5,2);
    
    if (r)
        % get red blob positions
        for i = 1:r
            [tempx,remain] = strtok(remain, ':*,|');
            [tempy,remain] = strtok(remain, ':*,|');
            red_data(i,1) = str2num(tempx);
            red_data(i,2) = str2num(tempy);
        end
    end
    if (g)
        % get green blob positions
        for i = 1:g
            [tempx,remain] = strtok(remain, ':*,|');
            [tempy,remain] = strtok(remain, ':*,|');
            green_data(i,1) = str2num(tempx);
            green_data(i,2) = str2num(tempy);
        end
    end
    if (b)
        % get blue blob positions
        for i = 1:b
            [tempx,remain] = strtok(remain, ':*,|');
            [tempy,remain] = strtok(remain, ':*,|');
            blue_data(i,1) = str2num(tempx);
            blue_data(i,2) = str2num(tempy);
        end    
    end
    
    [temp,remain] = strtok(remain, ':*,|#');

    time = str2num(temp);
    
    % store data in struct
    f1 = 'number_red';
    f2 = 'number_green';
    f3 = 'number_blue';
    v1 = r;
    v2 = g;
    v3 = b;
    
    f4 = 'red_data';
    f5 = 'green_data';
    f6 = 'blue_data';
    v4 = red_data;
    v5 = green_data;
    v6 = blue_data;
    
    f7 = 'time';
    v7 = time;
    
    f8 = 'cam_param';
    v8 = param;
    s = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7,f8,v8);
end
