cylinder87=[1.0642,1.0809,1.0809,1.0642,1.0157,1,0.9544,0.8169,0.6078,0.3065,0.2177,0.2246,0.2465,0.2543];

Re=[10000+4*9000:9000:172000-9000]



fur87=myData;
Re2=[10,15,20,25,30]*1.2045*D/1.815E-5


figure(1)
clf;
plot(Re2,fur87,'r','LineWidth',4);
hold on
plot(Re,cylinder87,'b','LineWidth',4);
legend('fur cylinder','smooth cylinder');

xlabel('Reynolds number');
ylabel('Drag coefficient');

