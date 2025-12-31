function img = RawFileLoaderMat(filename)

% this function needs comments
        % if the name "a3Tmp" change in future, this function will change
        % accordingly.
   
        display(['load ', filename]);
        
        a3Tmp = load(  filename );
        img = a3Tmp.a3ImageI16;
end