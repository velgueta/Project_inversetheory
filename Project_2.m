

%Velocity_up=rmoutliers(Data(:,11))
displa_up=[0.061 0.015 0.027 0.022 0.073 0.103]*0.01% displacement in m/yr
uncer_displa=[0.037 0.016 0.013 0.017 0.018 0.030]*0.01 % uncertantanties in m/yr
dd=displa_up';
mo=[50000000 30000000 5000000];%pa
mogicor=[475367.09 6580590.6];
%secondmogi=59.372121,-153.429222 coordenadas segunda morgi source
dosmogi=[475606.91 6581567.72];
tresmogi=[473435.58 6583934.51];
first=[479886.43 6581339.77];
second= [473801.83 6580062.95];
third=[475625.60 6577211.10];
fourth=[475126.14 6582593.14];
fifth=[474720.21 6580513.28];
sixth=[475971.04 6580583.54];

xest=[ 479886.43 473801.83 475625.60 475126.14 474720.21 475971.04 ];
xnorth=[6581339.77 6580062.95 6577211.10 6582593.14 6580513.28 6580583.54];
%deltap=model parameters
%calculate r, horizaontal distance
ro=zeros(length(xest),1) 
rd=zeros(length(xest),1)
rt=zeros(length(xest),1)
for h=1:length(xest)

    ro(h)=sqrt((xest(h)-mogicor(1))^2+(xnorth(h)-mogicor(2))^2)
     rd(h)=sqrt((xest(h)-dosmogi(1))^2+(xnorth(h)-dosmogi(2))^2)
      rt(h)=sqrt((xest(h)-tresmogi(1))^2+(xnorth(h)-tresmogi(2))^2)

end 




 %P= 3*10^6%$Pa$
 %a= 500 %$m$
  a= 200 %$m$

 GG= 30*10^9 %$Pa$
 d= 5000%m%
%%
 %Our model is gonna assume 3 differents sources of magma.
 
%Calculate G
   
%G= zeros(length(displa_up),1);
%for k=1:length(displa_up)
 %  G(k)=(3*d*(a^3))/(4*GG*(d^2+r(k)^2)^(3/2));
%end 
r=[ro rd rt];
G= zeros(length(displa_up),2);

for k=1:length(displa_up)
    for j=1:3
     G(k,j)=(3*d*(a^3))/(4*GG*(d^2+r(k,j)^2)^(3/2));
    end
end 
esquarevalues=[]
x=lsqnonneg(G,dd)
for y=1:length(d);
    dp=G*x;
    e=(d(y,1)-dp(y,1))
    esquare=e^2
    esquarevalues=[esquarevalues, esquare]
end
efinal=sum(esquarevalues)
l22=sqrt(efinal)

figure,
plot(ro,dd,'ko');
hold on
errorbar(r,dd,uncer_displa,'r*')
ylim([-0.003, 0.003])
legend('Dispacement rate + uncertainty ', 'FontSize', 12)
xlabel('Distance to the mogi source [m]',  'FontSize', 18)
ylabel('Displacement rate [m/year]', 'FontSize', 18)


figure,
subplot(3,3,1);
plot(ro,dd,'ko');
hold on
errorbar(r,dd,uncer_displa,'r*')
ylim([-0.003, 0.003])
legend('Dispacement rate + uncertainty ', 'FontSize', 3)
xlabel('Distance to the mogi source [m]',  'FontSize', 6)
ylabel('Displacement rate [m/year]', 'FontSize', 6)
subplot(3,3,2);
plot(rd,dd,'ko');
hold on
errorbar(r,dd,uncer_displa,'r*')
ylim([-0.003, 0.003])
legend('Dispacement + uncertainty ', 'FontSize', 3)
xlabel('Distance to the mogi source [m]',  'FontSize', 6)
ylabel('Displacement [m/year]', 'FontSize', 6)

subplot(3,3,3);
plot(rt,dd,'ko');
hold on
errorbar(r,dd,uncer_displa,'r*')
ylim([-0.003, 0.003])
legend('Dispacement + uncertainty ', 'FontSize', 3)
xlabel('Distance to the mogi source [m]',  'FontSize', 6)
ylabel('Displacement [m/year]', 'FontSize', 6)

%plot(rr,dd,'ko')

p=1;
[U S V] =svd(G)
Up = U(:,1:1);
Sp = S(1:1,1:1);
Vp = V(:,1:1);
% %(∆P) is pressure change
% %in the chamber (10-40 MPa),
  ginverse=Vp*inv(Sp)*Up'
  mest=ginverse*dd
% predict data using estimated model parameters
dpred = G*mest % [s] predicted data
resid = dd-dpred; % [s], residual between observed and predicted data
L2 = sqrt(resid'*resid)
R=Vp*Vp.'
unce_d_2=[uncer_displa(1,1)^2 uncer_displa(1,2)^2 uncer_displa(1,3)^2 uncer_displa(1,4)^2 uncer_displa(1,5)^2 uncer_displa(1,6)^2];%
disp("The units of unce_d_2 are (mm/yr)^2 ")
cov_d=diag(unce_d_2)
cov_m= ginverse*cov_d*(ginverse)'
uncertainities=[sqrt(cov_m(1,1)) sqrt(cov_m(2,2)) sqrt(cov_m(3,3))]
%%
disp('L2 calculate for the last row')

% plot data
figure
plot([dd dpred],'*');
legend('Observed data','Predicted data','FontSize', 15);
xlabel('Stations','FontSize', 18);
ylabel('Displacement rate [m/year]','FontSize', 18);
title(['Surface displacement up using S=1, L2 Norm = ',num2str(L2),' [m]']);

figure,imagesc(S)
title('Matrix S');
colorbar

figure,imagesc(R)
title('Matrix resolution using S=1');
colorbar

figure,imagesc(cov_m)
title('Covariance matrix using S=1');
colorbar
figure  % plot model parameters
subplot(2,1,1)
imagesc(reshape(mo,1,3)'); 
caxis([1 10]); axis equal off
colorbar
title('Idealized pressure [km/s]');
subplot(2,1,2)
imagesc(reshape(mest,1,3)'); 
title(['Change Presion, p=',num2str(p)]);
caxis([1 10]); axis equal off
colorbar