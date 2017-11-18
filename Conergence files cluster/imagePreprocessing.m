function [I0,I1] = imagePreprocessing(I0,I1)


sigma = 5;
I0 = imgaussfilt(I0,sigma);
I1 = imgaussfilt(I1,sigma);




end

