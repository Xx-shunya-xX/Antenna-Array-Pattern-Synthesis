function CF=SLL52(Amplitude)
theta1=[0:0.1:76];
theta2=[104:0.1:180];
theta=[theta1,theta2];
bb=1:length(theta);
CF=max(20*log10(abs(2*Amplitude(1) * cos(pi/2*cos(pi/180*theta(bb))) + 2*Amplitude(2) * cos(3*pi/2*cos(pi/180*theta(bb))) + 2*Amplitude(3)*cos(5*pi/2*cos(pi/180*theta(bb))) + 2*Amplitude(4) * cos(7*pi/2*cos(pi/180*theta(bb))) + 2*Amplitude(5) * cos(9*pi/2*cos(pi/180*theta(bb))) )));
end
