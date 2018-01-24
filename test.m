f = fopen('salida_real', 'r');
%g = fopen('salida_imaginaria', 'r');
R = fread(f, 'double');
%I = fread(g, 'double');
fclose(f);
%fclose(g);
R = reshape(R,512,512);
%I = reshape(I,512,512);
%H = R + I*1j;
%H = H';
H = R';

V = ifft2(H);
F = ifftshift(V);
imagesc(abs(F));axis('square');colormap(jet);
figure
imagesc(log(abs(V)+1));axis('square');colormap(gray);
        