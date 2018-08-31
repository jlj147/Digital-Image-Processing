%Jayce Jones
%EE 4323 Spring 2018
%Assignment_3

function fftshow(f,type)

% Displays the fft matrix F using imshow, where TYPE must be one of
% 'abs' or 'log'. If TYPE='abs', then then abs(f) is displayed; if
% TYPE='log' then log(1+abs(f)) is displayed. If TYPE is omitted, then
% 'log' is chosen as a default.

if nargin<2
type='log';
end

if (type=='log')
fl = log(1+abs(f));
fm = max(fl(:));
imshow(im2uint8(fl/fm))

elseif (type=='abs')
fa=abs(f);
fm=max(fa(:));
imshow(fa/fm)
else
error('TYPE must be abs or log.');
end



