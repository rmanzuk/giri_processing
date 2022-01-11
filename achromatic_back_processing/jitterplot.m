function [] = jitterplot(values,labels,jitter,colors,groups)
    %% begin the function
    markers = {'o','s','d','^','p'};
    % we'll just loop through the number of classes and plot everything
    hold on
    for i = 1:max(labels)
        % good to know the number of points in this class
        n_points = sum(labels == i);
        
        % make some random values that will account for the jitter
        for_jitter = rand(1,n_points);
        
        % scale those values to match the level of jitter desired
        for_jitter = (for_jitter*2*jitter)-jitter;
        
        % now we're ready to scatter x is just i plus jitter and y are the
        % values
        these_values = values(labels==i);
        
        scatter(i + for_jitter,these_values,30,colors(i,:),'filled',markers{groups(i)},'MarkerEdgeColor',[0,0,0])
    end
    hold off
end