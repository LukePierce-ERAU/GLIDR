free = struct();
pull = struct();
glide = struct();
steady = struct();
master = struct();


for free_ind = 1:size(freevar,3)
    free_mat_ind = cell2mat(freevar(:,:,free_ind));
    t = free_mat_ind(1:size(free_mat_ind,1)/3);
    x = zeros(size(free_mat_ind,1)/3,1);
    alt = free_mat_ind(size(free_mat_ind,1)/3+1:size(free_mat_ind,1)/3*2);
    speed = free_mat_ind(size(free_mat_ind,1)/3*2+1:end);
    t_name = ['t', num2str(free_ind)];
    x_name = ['x', num2str(free_ind)];
    alt_name = ['alt', num2str(free_ind)];
    speed_name = ['speed', num2str(free_ind)];
    free.(t_name) = t;
    free.(x_name) = x;
    free.(alt_name) = alt;
    free.(speed_name) = speed;
end

for pull_ind = 1:size(pullvar,3)
    pull_mat_ind = cell2mat(pullvar(:,:,pull_ind));
    t = free_out(1,1,pull_ind) + pull_mat_ind(1:size(pull_mat_ind,1)/5);
    x = pull_mat_ind(size(pull_mat_ind,1)/5+1:size(pull_mat_ind,1)/5*2);
    alt = pull_mat_ind(size(pull_mat_ind,1)/5*2+1:size(pull_mat_ind,1)/5*3);
    speed = pull_mat_ind(size(pull_mat_ind,1)/5*3+1:size(pull_mat_ind,1)/5*4);
    t_name = ['t', num2str(pull_ind)];
    x_name = ['x', num2str(pull_ind)];
    alt_name = ['alt', num2str(pull_ind)];
    speed_name = ['speed', num2str(pull_ind)];
    pull.(t_name) = t;
    pull.(x_name) = x;
    pull.(alt_name) = alt;
    pull.(speed_name) = speed;

end

for glide_ind = 1:size(glidevar,3)
    glide_mat_ind = cell2mat(glidevar(:,:,glide_ind));
    t = free_out(1,1,glide_ind) + pull_out(1,1,glide_ind) + glide_mat_ind(1:size(glide_mat_ind,1)/4);
    x = glide_mat_ind(size(glide_mat_ind,1)/4+1:size(glide_mat_ind,1)/4*2);
    alt = glide_mat_ind(size(glide_mat_ind,1)/4*2+1:size(glide_mat_ind,1)/4*3);
    speed = glide_mat_ind(size(glide_mat_ind,1)/4*3+1:end);
    t_name = ['t', num2str(glide_ind)];
    x_name = ['x', num2str(glide_ind)];
    alt_name = ['alt', num2str(glide_ind)];
    speed_name = ['speed', num2str(glide_ind)];
    glide.(t_name) = t;
    glide.(x_name) = x;
    glide.(alt_name) = alt;
    glide.(speed_name) = speed;

end

for steady_ind = 1:size(steadyvar,3)
    steady_mat_ind = cell2mat(steadyvar(:,:,steady_ind));
    t = steady_inip(1,steady_ind) + steady_mat_ind(1:size(steady_mat_ind,1)/4);
    x = steady_mat_ind(size(steady_mat_ind,1)/4+1:size(steady_mat_ind,1)/4*2);
    alt = steady_mat_ind(size(steady_mat_ind,1)/4*2+1:size(steady_mat_ind,1)/4*3);
    speed = steady_mat_ind(size(steady_mat_ind,1)/4*3+1:end);
    t_name = ['t', num2str(steady_ind)];
    x_name = ['x', num2str(steady_ind)];
    alt_name = ['alt', num2str(steady_ind)];
    speed_name = ['speed', num2str(steady_ind)];
    steady.(t_name) = t;
    steady.(x_name) = x;
    steady.(alt_name) = alt;
    steady.(speed_name) = speed;
end

% Creating one line for each config

for jj = 1:1:ii
master.(['t', num2str(jj)]) = vertcat(free.(['t',num2str(jj)]),pull.(['t',num2str(jj)]),glide.(['t',num2str(jj)]),steady.(['t',num2str(jj)]));
master.(['x', num2str(jj)]) = vertcat(free.(['x',num2str(jj)]),pull.(['x',num2str(jj)]),glide.(['x',num2str(jj)]),steady.(['x',num2str(jj)]));
master.(['alt', num2str(jj)]) = vertcat(free.(['alt',num2str(jj)]),pull.(['alt',num2str(jj)]),glide.(['alt',num2str(jj)]),steady.(['alt',num2str(jj)]));
master.(['speed', num2str(jj)]) = vertcat(free.(['speed',num2str(jj)]),pull.(['speed',num2str(jj)]),glide.(['speed',num2str(jj)]),steady.(['speed',num2str(jj)]));
end

figure
hold on
plot(master.t1/60,master.alt1/1000,'g')
plot(master.t2/60,master.alt2/1000,'k')
plot(master.t3/60,master.alt3/1000,'b')
xline(90,'r')
grid on
title('Altitude over Time')
xlabel('Time [mins]')
ylabel('Altitude [km]')

figure
hold on
plot(master.x1/1000,master.alt1/1000,'g')
plot(master.x2/1000,master.alt2/1000,'k')
plot(master.x3/1000,master.alt3/1000,'b')
xline(45,'r')
grid on
title('2D Flight path')
xlabel('Distance Travelled [km]')
ylabel('Altitude [km]')

figure
hold on
plot(master.speed1,master.alt1/1000,'g')
plot(master.speed2,master.alt2/1000,'k')
plot(master.speed3,master.alt3/1000,'b')
grid on
title('Altitude over Speed')
xlabel('Speed [m/s]')
ylabel('Altitude [km]')