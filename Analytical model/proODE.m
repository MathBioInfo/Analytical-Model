function dy = proODE(t,y)
dy = zeros(8,1);
global rL rS rD rI;

%rL is relative rate of Lysogeny
%rS is the relative rate of selective benfit
%rD is the relative rate of degraddation
%rI is the relative rate of induction

%Definition of the ODE system

A = (rL + rS - rI).*y(1) + (rL - rI).*y(2) + (rS).*y(3) + (rS - rI).*y(4) + (-rI).*y(6) + (rS).*y(7);

dy(1) = (rL + rS -3*rD - rI )*y(1) - A*y(1); %y(1) contain all genes
dy(2) = (rL-2*rD-rI)*y(2)+rD*y(1) - A*y(2);   %y(2) contain lys and inf genes and no ben genes
dy(3) = (rS-2*rD)*y(3)+rD*y(1) - A*y(3);      %y(3) contain ben and inf genes and no lys genes
dy(4) = (rS-2*rD-rI)*y(4)+rD*y(1) - A*y(4);   %y(4) contain ben and lys genes and no inf genes
dy(5) = (-rD)*y(5)+rD*y(2)+rD*y(3) - A*y(5);  %y(5) contain only inf genes and no ben, lys genes
dy(6) = (-rD-rI)*y(6)+rD*y(2)+rD*y(4) - A*y(6); %y(6) contain only lys genes and no ben, ing genes
dy(7) = (rS-rD)*y(7)+rD*y(3)+rD*y(4) - A*y(7);  %y(7) contain only ben genes and no lys or inf genes
dy(8) = rD*(y(5) + y(6) + y(7)) - A*y(8);
 %%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%
