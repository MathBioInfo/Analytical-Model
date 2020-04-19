% MATLAB code for solving prophage simulation model
clear, clear all;
global rL rS rD rI;
GeneType = 3;
influx = 0.0;
rS = 0.52;
rL = 1.2;
rD = 0.01;
rI = 1;
y0 = [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]'; % inatial conditions for solving the system
options = odeset('RelTol', 1e-4, 'Abs', 1e-8*[1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7]);
[T,Y] = ode45(@proODE, [0 200], y0, options);
negi = length(T)*ones(1,8);
for i=1:8
  negs = find(Y(:,i)<0);
  if ~isempty(negs) negi(i) = min(negs); end
end
imax = min(negi)-1;
fprintf(1,'Numerical error: a population was negative at time t=%f\n',T(imax+1));
T = T(1:imax);
Y = Y(1:imax,:);
pro(:,1) = Y(:,1);
pro(:,2) = Y(:,2);
pro(:,3) = Y(:,3);
pro(:,4) = Y(:,4);
pro(:,5) = Y(:,5);
pro(:,6) = Y(:,6);
pro(:,7) = Y(:,7);
pro(:,8) = Y(:,8);
prosum(:,1) = 1.*pro(:,1)+ 1.*pro(:,2)+1.*pro(:,3)+ 1.*pro(:,4)+1.*pro(:,5)+1.*pro(:,6)+1.*pro(:,7) + pro(:,8);
ave = zeros(length(T), GeneType);
ave(:,1) = (pro(:,1) + pro(:,3) + pro(:,4) + pro(:,7))./(prosum(:,1)); % beneficial genes
ave(:,2) = (pro(:,1) + pro(:,2) + pro(:,4) + pro(:,6))./(prosum(:,1)); % lysis genes
ave(:,3) = (pro(:,1) + pro(:,2) + pro(:,3) + pro(:,5))./(prosum(:,1)); % inf genes
figure(2)
plot(T, ave(:,1), 'g', 'LineWidth',2); %ben genes
hold on;
plot(T, ave(:,2), 'r', 'LineWidth',2); %lys genes
hold on;
plot(T, ave(:,3), 'y', 'LineWidth',2); %inf genes
legend({'Beneficial genes', 'Lysis genes', 'Infectious genes'},'Interpreter','latex');
xlabel({'Time (Prophage life time)'},'Interpreter','latex')
ylabel({'Average number of genes of each type'},'Interpreter','latex')
hold off;
