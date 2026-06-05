figure

plot(rad2deg(RQA.theta),RQA.rOuterN,'.')
hold on

plot(rad2deg(RQA.theta),RQA.rOuterFit,'LineWidth',2)

xline(rad2deg(theta_bad),'r')

%% ...
figure
plot(rad2deg(RQA.theta),RQA.nYSegments,'.-')
ylabel('Number of Segments')