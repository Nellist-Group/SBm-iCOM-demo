function sel_range = RangeFunc(img)
% this function needs comments

   figure;
   [~,RECT] = imcrop(img./max(img(:)));
   sel_range = round(RECT);
   close; 
end