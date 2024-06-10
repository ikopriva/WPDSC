function swp2show(swpc)

lsum = size(swpc,3);

% level 1
c = 'AHVD';
figure,
for i=1:4,
    s = sprintf('%d) %c', i, c(i));
    subplot(2,2,i), imshow(swpc(:,:,i), []), title(s)
end


% level 2
ind = [1 2 5 6 3 4 7 8 9 10 13 14 11 12 15 16];
c = 'AAAHAVADHAHHHVHDVAVHVVVDDADHDVDD';
k = 5;
if k  > lsum return, end
figure,
for i = 1:16 
    if (k+i-1) > lsum break, end
    s = sprintf('%d) %s', i, c(2*(i-1)+[1 2]));
    subplot(4,4,ind(i)), imshow(swpc(:,:,k+i-1), []), title(s)
end

% the rest
for k = 21:16:lsum 
    if k  > lsum return, end
    figure,
    for i = 1:16
        if (k+i-1) > lsum break, end
        subplot(4,4,ind(i)), imshow(swpc(:,:,k+i-1), []), title(int2str(k+i-1))
    end
end
