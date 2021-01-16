% The method is descried in detail in Main text and Supplementary
% information. We have opted not to make all computer code publicly
% available, because the code will be used for other projects. We would
% like to provide the procedure to process data.
% If further information needed, please contact us.
% 01/16/2021  for  publication in Nature Communication
clear all;
close all;

data = dmread('EELS Spectrum Image_MOF.dm4');   % read dm4 data
 
I = data;          % load image 
Im = data;                   
for i =1:1:size(I,1); 
    for j=1:size(I,2);
    Im(i,j,:) = Im(i,j,:) - mean(Im(i,j,1:20));
    Im1(i,j) = sum(Im(i,j,:));
    end
end 

xaxis = -9.35:0.05:93;

Data = normZLP(data,xaxis,1,true);    % Aligh ZLP and normalize spectra
I1 = Data;
for i =1:1:size(I1,1); 
    for j=1:size(I1,2);
    imm(i,j) = sum(I1(i,j,:));      
    end
end 

S(:)= Data(9,25,:);   % plot the spectrum 
figure,
line1 = plot(xaxis,S,'k','linewidth',2);
hold on

II = load('compandim_1to65_comp3.mat'); % load NMF results
comp1 = II.wproaglspcr(:,1);  % first component
C1 = comp1/max(comp1);
s1 = plot(II.xaxisproaglspcr, C1,'b','linewidth',2);    % plot first comp1
hold on

% fitting like gaussian fitting, double gaussianfitting can be done here


% format of plot
ylabel('\fontsize{16}Normalized intensity');
xlabel('\fontsize{16}Energy loss (eV)');
set(gca,'tickdir','out'); 
xlim([1 50]);
set(gca,'fontsize',15);
set(gca,'XTick',0:5:50)
ylim([-0.0005 1.2]);

xL=xlim;yL=ylim;
plot(xL,[yL(2),yL(2)],'k',[xL(2),xL(2)],[yL(1),yL(2)],'k','LineWidth',1.2)
box off
axis([xL yL]) ;
set(gca,'linewidth',1.2);

yticks([''])

l = legend([s1,line1],'MOF(l)','raw data');
set(l,'box','off');
set(l, 'FontSize', 14);
