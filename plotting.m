close all;
figure(1)
t2 = [tseg fliplr(tseg)]; % Prepare for plotting bounds via "fill"
x2 = [xhatseg-boundseg fliplr(xhatseg+boundseg)];
h1a = fill(t2,x2(1,:),'b',"FaceAlpha",0.05); hold on; grid on;
fill(t2,x2(2,:),'b',"FaceAlpha",0.05);

h2 = plot(tseg,xseg',tseg,xhatseg','--'); ylim([-0.13 0.13]);
legend([h2;h1a],{'True position','True velocity','AEKF Position estimate',...
  'AEKF Velocity estimate','Confidence bounds'},'Location','BestOutside');
title('Demonstrating state estimates for AEKF');
xlabel('Time (s)'); ylabel('State (m or m/s)');

figure(2)
fill([tseg fliplr(tseg)],[-boundseg(1,:) fliplr(boundseg(1,:))],'b',...
  "FaceAlpha",0.05); hold on; grid on;
plot(tseg,errseg(1,:),'b');
title('Position estimation error for AEKF');
xlabel('Time (s)'); ylabel('Error (m)');

figure(3)
fill([tseg fliplr(tseg)],[-boundseg(2,:) fliplr(boundseg(2,:))],'b',...
  "FaceAlpha",0.05); hold on; grid on;
plot(tseg,errseg(2,:),'b'); ylim([-0.05 0.05]);
title('Velocity estimation error for AEKF');
xlabel('Time (s)'); ylabel('Error (m/s)');

% I plot every tenth point for speed
figure(4)
h1 = plot(t(1:10:end)/3600,Qstore(1,1:10:end)); title('Evolution of estimate of SigmaW'); grid on
xlabel('Time (h)');
ylabel('SigmaW values'); hold on
h2 = plot(t(1:10:end)/3600,Qstore(3:4,1:10:end));
xlim([0 3]);
legend([h1;h2],{'(1,1) component','(1,2)=(2,1) component','(2,2) component'},'location','bestOutside');


figure(5)
h1=plot(t(1:10:end)/3600,Rstore(1:10:end)); title('Evolution of estimate of SigmaV'); grid on
xlabel('Time (h)'); ylabel('SigmaV values')
ylim([0 1.3e-3]); xlim([0 3]);
