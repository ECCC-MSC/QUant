function drawBoxVert(draw_data, lineWidth, xpos, width, colour, logScale)

n = size(draw_data, 2);

% for log scale plots, need a width that depends on the values of the box

unit = width/2;


for i = 1:n
    
    if logScale
        unit = xpos(i)*.1;
    end
    
    v = draw_data(:,i);
    
    % draw the min line
    plot([xpos(i)-unit, xpos(i)+unit], [v(5), v(5)], 'LineWidth', lineWidth, 'Color', colour);
    % draw the max line
    plot([xpos(i)-unit, xpos(i)+unit], [v(1), v(1)], 'LineWidth', lineWidth, 'Color', colour);
    % draw middle line
    plot([xpos(i)-unit, xpos(i)+unit], [v(3), v(3)], 'LineWidth', lineWidth, 'Color', colour);
    % draw vertical line
    plot([xpos(i), xpos(i)], [v(5), v(4)], 'LineWidth', lineWidth, 'Color', colour);
    plot([xpos(i), xpos(i)], [v(2), v(1)], 'LineWidth', lineWidth, 'Color', colour);
    % draw box
    plot([xpos(i)-unit, xpos(i)+unit, xpos(i)+unit, xpos(i)-unit, xpos(i)-unit], [v(2), v(2), v(4), v(4), v(2)], 'LineWidth', lineWidth, 'Color', colour);
end