function [t1,t2,t3,t4] = getlines(fd1,fd2,fd3,fd4)
    t1 = fgetl(fd1);
    t2 = fgetl(fd2);
    t3 = fgetl(fd3);
    t4 = fgetl(fd4);
end