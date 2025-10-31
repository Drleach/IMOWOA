function PlotCosts(pop)

    Costs=[pop.Cost];
    
    plot(Costs(1,:),Costs(2,:),'r*','MarkerSize',8);
    xlabel('Objective1');
    ylabel('Objective2');
    title('Non-dominated solution');
    grid on;
    set(gcf,'color','w')
end