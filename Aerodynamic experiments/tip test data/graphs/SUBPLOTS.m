%CD v Reynolds

h1=openfig('CDvR0.fig','reuse');
ax1=gca;
h2=openfig('CDvR15.fig','reuse');
ax2=gca;
h3=openfig('CDvR30.fig','reuse');
ax3=gca;
h4=openfig('CDvR45.fig','reuse');
ax4=gca;
h5=openfig('CDvR60.fig','reuse');
ax5=gca;

h6=figure;


s1=subplot(3,2,1,'XScale','log');
title('Cd vs Reynolds number at 0 degrees')
xlabel('Reynolds number')
ylabel('Cd')
s2=subplot(3,2,2,'XScale','log');
title('Cd vs Reynolds number at 15 degrees')
xlabel('Reynolds number')
ylabel('Cd')
s3=subplot(3,2,3,'XScale','log');
title('Cd vs Reynolds number at 30 degrees')
xlabel('Reynolds number')
ylabel('Cd')
s4=subplot(3,2,4,'XScale','log');
title('Cd vs Reynolds number at 45 degrees')
xlabel('Reynolds number')
ylabel('Cd')
s5=subplot(3,2,[5,6],'XScale','log');
title('Cd vs Reynolds number at 60 degrees')
xlabel('Reynolds number')
ylabel('Cd')


fig1=get(ax1,'children');
fig2=get(ax2,'children');
fig3=get(ax3,'children');
fig4=get(ax4,'children');
fig5=get(ax5,'children');

suptitle('Cd vs Reynolds number at different angles of attack')

copyobj(fig1,s1); 
copyobj(fig2,s2);
copyobj(fig3,s3);
copyobj(fig4,s4);
copyobj(fig5,s5);
legend('fur','spray','cheetah','cheetah spray');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Cd v angle
h7=openfig('CDvAbf.fig','reuse');
ax7=gca;
h8=openfig('CDvAbfs.fig','reuse');
ax8=gca;
h9=openfig('CDvAc.fig','reuse');
ax9=gca;
h10=openfig('CDvAcs.fig','reuse');
ax10=gca;

h11=figure;

s7=subplot(2,2,1);
title('Cd vs Angle for black fur')
xlabel('Angle (degrees)')
ylabel('Cd')
s8=subplot(2,2,2);
title('Cd vs Angle for black fur with spray')
xlabel('Angle (degrees)')
ylabel('Cd')
s9=subplot(2,2,3);
title('Cd vs Angle for cheetah fur')
xlabel('Angle (degrees)')
ylabel('Cd')
s10=subplot(2,2,4);
title('Cd vs Angle for cheetah fur with spray')
xlabel('Angle (degrees)')
ylabel('Cd')


fig7=get(ax7,'children');
fig8=get(ax8,'children');
fig9=get(ax9,'children');
fig10=get(ax10,'children');

suptitle('Cd vs Angle of attack for various speeds')

copyobj(fig7,s7); 
copyobj(fig8,s8);
copyobj(fig9,s9);
copyobj(fig10,s10);
legend('10m/s','15m/s','20m/s','25m/s','30m/s');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %CD v speed
h12=openfig('CDvSbf.fig','reuse');
ax12=gca;
h13=openfig('CDvSbfs.fig','reuse');
ax13=gca;
h14=openfig('CDvSc.fig','reuse');
ax14=gca;
h15=openfig('CDvScs.fig','reuse');
ax15=gca;

h16=figure;

s12=subplot(2,2,1);
title('Cd vs Speed for black fur')
xlabel('Speed (m/s)')
ylabel('Cd')
s13=subplot(2,2,2);
title('Cd vs Speed for black fur with spray')
xlabel('Speed (m/s)')
ylabel('Cd')
s14=subplot(2,2,3);
title('Cd vs Speed for cheetah fur')
xlabel('Speed (m/s)')
ylabel('Cd')
s15=subplot(2,2,4);
title('Cd vs Speed for cheetah fur with spray')
xlabel('Speed (m/s)')
ylabel('Cd')


fig12=get(ax12,'children');
fig13=get(ax13,'children');
fig14=get(ax14,'children');
fig15=get(ax15,'children');

suptitle('Cd vs speed at different angles of attack')

copyobj(fig12,s12); 
copyobj(fig13,s13);
copyobj(fig14,s14);
copyobj(fig15,s15);
legend('0','15','30','45','60');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%force plots 
h17=openfig('FvSbf.fig','reuse');
ax17=gca;
h18=openfig('FvSbfs.fig','reuse');
ax18=gca;
h19=openfig('FvSc.fig','reuse');
ax19=gca;
h20=openfig('FvScs.fig','reuse');
ax20=gca;

h21=figure;

s17=subplot(2,2,1);
title('F vs Speed for black fur')
xlabel('Speed (m/s)')
ylabel('F (N)')
s18=subplot(2,2,2);
title('F vs Speed for black fur with spray')
xlabel('Speed (m/s)')
ylabel('F (N)')
s19=subplot(2,2,3);
title('F vs Speed for cheetah fur')
xlabel('Speed (m/s)')
ylabel('F (N)')
s20=subplot(2,2,4);
title('F vs Speed for cheetah fur with spray')
xlabel('Speed (m/s)')
ylabel('F (N)')


fig17=get(ax17,'children');
fig18=get(ax18,'children');
fig19=get(ax19,'children');
fig20=get(ax20,'children');

suptitle('Force vs Speed at different angles of attack')

copyobj(fig17,s17); 
copyobj(fig18,s18);
copyobj(fig19,s19);
copyobj(fig20,s20);
legend('0','15','30','45','60');