%Calculation of relative errors Cd
close all
open("matlab.mat")
%Calculation of the relative error of CFD with respect to Missile Datcom
%Case Rough Surface
Rel_Rough = abs((CDcfdfinal(:,2)-cdnAIRROUGH(:,2)))./cdnAIRROUGH(:,2);
Rel_Smooth= abs((CDcfdfinal(:,2)-cdnAIRSMOOTH(:,2)))./cdnAIRSMOOTH(:,2);


figure(1)
Rel_Rough = abs((CDcfdfinal(:,2)-cdnAIRROUGH(:,2)))./cdnAIRROUGH(:,2);
plot(CDcfdfinal(:,1),Rel_Rough,'r', 'LineWidth',2)
hold on 
plot(CDcfdfinal(:,1),CDcfdfinal(:,2),'b', 'LineWidth',2)
hold on
plot(cdnAIRROUGH(:,1),cdnAIRROUGH(:,2),'k--', 'LineWidth',2)
hold on
plot(cdnAIRSMOOTH(:,1),cdnAIRSMOOTH(:,2),'m--', 'LineWidth',2)
legend("Relative error of cfd with respect to missile datcom simulation with rough surface", "Cd CFD simulations", "Cd Datcom Simulations with rough surface","Cd Datcom Simulations with smooth surface")


figure(2)
Rel_Smooth= abs((CDcfdfinal(:,2)-cdnAIRSMOOTH(:,2)))./cdnAIRSMOOTH(:,2);
plot(CDcfdfinal(:,1),Rel_Smooth,'r', 'LineWidth',2)
hold on 
plot(CDcfdfinal(:,1),CDcfdfinal(:,2),'b', 'LineWidth',2)
hold on
plot(cdnAIRROUGH(:,1),cdnAIRROUGH(:,2),'k--', 'LineWidth',2)
hold on
plot(cdnAIRSMOOTH(:,1),cdnAIRSMOOTH(:,2),'m--', 'LineWidth',2)
legend("Relative error of cfd with respect to missile datcom simulation with smooth surface", "Cd CFD simulations", "Cd Datcom Simulations with rough surface","Cd Datcom Simulations with smooth surface")

figure(3)
Rel_Smooth= abs((CDcfdfinal(:,2)-cdnAIRSMOOTH(:,2)))./cdnAIRSMOOTH(:,2);
plot(CDcfdfinal(:,1),Rel_Smooth,'r', 'LineWidth',2)
hold on
Rel_Rough = abs((CDcfdfinal(:,2)-cdnAIRROUGH(:,2)))./cdnAIRROUGH(:,2);
plot(CDcfdfinal(:,1),Rel_Rough,'b', 'LineWidth',2)
legend("Relative error of CFD Simulation with respect to MD simulation with smooth surface","Relative error of CFD Simulation with respect to MD simulation with rough surface")
