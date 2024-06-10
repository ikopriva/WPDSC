function showFace(y)
    close all;
%    subplot(1,2,1);
    imshow(reshape(y,48,42),[min(y) max(y)]);
%    subplot(1,2,2);
%    imshow(reshape(-y,48,42),[-max(y) -min(y)]);
end