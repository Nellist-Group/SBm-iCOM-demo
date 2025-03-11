function [ output] = RoadmapScanner( RonchiRoadmap, SelImFrame, sel_range  )
% this function needs comments


   px = [sel_range(1), sel_range(1)+sel_range(3)-1];
   py = [sel_range(2), sel_range(2)+sel_range(4)-1];
   
for i=1:length(RonchiRoadmap.Flag)
    
    if RonchiRoadmap.Flag(i) == SelImFrame
        
        if RonchiRoadmap.Xcoord(i) > px(2)&& RonchiRoadmap.Ycoord(i) > py(2)
            break;
        elseif RonchiRoadmap.Xcoord(i) <= px(2)&& RonchiRoadmap.Ycoord(i) <= py(2)...
                && RonchiRoadmap.Xcoord(i) >= px(1) && RonchiRoadmap.Ycoord(i) >= py(1)
            
            if (exist('fname' , 'var') == 0) || ~strcmp(RonchiRoadmap.Filename{i},fname)
                fname = RonchiRoadmap.Filename{i};
                img = RawFileLoaderMat( fname );
            end
                       
            % This is where some task happens:
            xx = RonchiRoadmap.Xcoord(i)-px(1)+1;
            yy = RonchiRoadmap.Ycoord(i)-py(1)+1;
            output(:,:,yy,xx) = img(:,:,RonchiRoadmap.Seq(i));
            %display([num2str(RonchiRoadmap.Xcoord(i)),'/',num2str(RonchiRoadmap.Ycoord(i))])
            
        end
    end
end

if ~isa(class(output),'double')
    output = double(output);
end

end

