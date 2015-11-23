%% Visualises the paths taken by the NN-powered Sheepdog

clear

stem = 'log';

i = 1;
while 1
    if i < 10
        filenames(i,:) = strcat(stem,'0',num2str(i),'.mat');
    else
        filenames(i,:) = strcat(stem,num2str(i),'.mat');
    end
    if exist(filenames(i,:),'file')
        i = i + 1;
    else
        filenames(i,:) = [];
        break
    end
end

figure
axis([0 500 0 500]);
title(stem);
hold on;
for i = 1:size(filenames,1)
    load(filenames(i,:));
    plot(history.mouse_pos(:,1),history.mouse_pos(:,2));
end
hold off;