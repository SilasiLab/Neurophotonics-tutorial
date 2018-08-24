function  ImageROIMean = MeanFromROI(Image,XUpperLeft,YUpperLeft,ROISize)

%X and Y inputs must be reversed and adjusted by 1 elsewhere
roi    = Image(YUpperLeft:YUpperLeft + ROISize-1,XUpperLeft:XUpperLeft + ROISize-1,:);
signal = mean(mean(roi));
signal = double(signal);
ImageROIMean = reshape(signal,1,length(signal));

end