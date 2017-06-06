function visualize_ced(dist, options)
    nData = length(dist);
    x = 0:0.05:1;
    x = x*(max(dist) - min(dist)) + min(dist);
    CED = zeros(length(x),1);
    c = 0;

    for thres = x
        c = c + 1;
        idx = find(dist <= thres);
        CED(c) = length(idx)/nData;
    end

    
    plot( x, CED, 'LineWidth', 2 , 'MarkerEdgeColor','r');
    axis([min(dist) max(dist) 0 1])
    title(['CED: Dist average: ' num2str(mean(dist))]);
    grid on;
    saveas(gcf,[options.ResultFigurePath 'CED.jpg']);
end