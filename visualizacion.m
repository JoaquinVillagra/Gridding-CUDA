f = fopen('HLTauB6.raw', 'r');
I = fread(f, 'double');
fclose(f);
F = reshape(I,12785504/4,4);
%I = I';
%V = fft2(I);
%F = fftshift(V);
%imagesc(I);axis('square');colormap(gray);
%figure
%imagesc(log(abs(V)+1));axis('square');colormap(gray);
%e=I(:,[1,2]);
k = I(:,[3,4]);
%scatter(e(:,1),e(:,2));
%e = e * 0.00000484814;
