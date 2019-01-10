function [ out ] = ROI2indx( in )
%ROI2rect Takes ROI type coordinates (xmin ymin width height) and converts them to index  type coordinate
%(rowmin,rowmax,colmin,colmax)

out(1)=in(2);
out(2)=in(2)+in(4);
out(3)=in(1);
out(4)=in(1)+in(3);

end

