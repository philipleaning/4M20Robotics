dim = size(history.sheep_x,1);

history.grid_status = zeros(dim,25);

for i = 1:dim
    for y = 100:100:500
        for x = 100:100:500
            for j = 1:size(history.sheep_x,2)
                    if (history.sheep_y(i,j) > y-100) && (history.sheep_y(i,j) < y) && (history.sheep_x(i,j) > x-100) && (history.sheep_x(i,j) < x)
                        history.grid_status(i,(x/100)+(y/20 - 5)) = 1;
                    end
            end
        end
    end
end

clearvars -except history