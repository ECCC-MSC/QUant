function drawBoxVert(draw_data, lineWidth, ypos, width, color, logScale)

n = size(draw_data, 2);

unit = width/2;
    
for i = 1:n
    
        
    if logScale
        unit = ypos(i)*.05;
    end
    
    v = draw_data(:,i);
    
    % draw the min line
    plot([v(5), v(5)], [ypos(i)-unit, ypos(i)+unit],  'LineWidth', lineWidth,  'Color', color);
    % draw the max line
    plot([v(1), v(1)], [ypos(i)-unit, ypos(i)+unit], 'LineWidth', lineWidth,  'Color', color);
    % draw middle line
    plot([v(3), v(3)], [ypos(i)-unit, ypos(i)+unit], 'LineWidth', lineWidth,  'Color', color);
    % draw vertical line
    plot([v(5), v(4)],[ypos(i), ypos(i)],  'LineWidth', lineWidth,  'Color', color);
    plot([v(2), v(1)],[ypos(i), ypos(i)],  'LineWidth', lineWidth,  'Color', color);
    % draw box
    plot([v(2), v(2), v(4), v(4), v(2)], [ypos(i)-unit, ypos(i)+unit, ypos(i)+unit, ypos(i)-unit, ypos(i)-unit], 'LineWidth', lineWidth, 'Color', color);
end