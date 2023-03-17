function den = ves_density(img)
%   calculate vessel density
%   Input: binary vessel segments
%   Jiarui Yang
img(img(:)~=0)=1;
den = sum(img(:))/numel(img);
end

