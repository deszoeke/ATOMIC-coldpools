function inset_plot(i,tcold,Taf,cp_matrix,qairf,dDf,rrf)
    axes('Position',[.69 .55 .15 .15])
    box on
    plot(tcold(i,1:61),Taf(cp_matrix(i,1:61)),'.');
    datetick('x','HH:MM','keeplimits','keepticks')
    ylabel('Ta [C]')
    title(['CP #',num2str(i),'; onset on ',datestr(tcold(i,1))])
%     ylim([22.5 28])

    axes('Position',[.69 .35 .15 .15])
    box on
    plot(tcold(i,1:61),qairf(cp_matrix(i,1:61)),'.');
    datetick('x','HH:MM','keeplimits','keepticks')
    ylabel('q [g/kg]')
%     ylim([11 18])

    axes('Position',[.69 .15 .15 .15])
    box on
    plot(tcold(i,1:61),dDf(cp_matrix(i,1:61)),'.');
    datetick('x','HH:MM','keeplimits','keepticks')
    ylabel('\deltaD [permil]')
    yyaxis right
    plot(tcold(i,1:61),rrf(cp_matrix(i,1:61)),'.r')
    ylabel('RR [mm/hr]')

%     ylim([-80 -64])