f = fopen('HLreduced.raw', 'r');
I = fread(f, 'float');
fclose(f);
I = reshape(I,16,16);
I = I';
V = fft2(I);
F = fftshift(fftshift(V));
imagesc(I);axis('square');colormap(gray);
figure
imagesc(log(abs(V)+1));axis('square');colormap(gray);