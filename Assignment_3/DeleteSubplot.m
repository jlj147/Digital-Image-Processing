%Jayce Jones
%EE 4323 Spring 2018
%Assignment_3

function DeleteSubplot()

%Deletes subplots that have been created by other enhancements when the
%subplots no longer need to be dispalyed

SubA=subplot(2,2,3);
SubB=subplot(2,2,4);
delete(SubA);
delete(SubB);

end

