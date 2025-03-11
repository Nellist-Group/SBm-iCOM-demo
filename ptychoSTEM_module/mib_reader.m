function C = ReadImage(Filename);
% Returns Merlin image from .mib file
fid =fopen(Filename,'r','b');
HDR = char(fread(fid,50));      % get first 50 bytes of header
firstHDR = textscan(HDR, '%s %d %d %d %d %d  %s %s','delimiter', ',');  % parse into values
s = cell2mat(firstHDR(3));        % data offset
chips = cell2mat(firstHDR(4));  % Number of chips
x = cell2mat(firstHDR(5));      % image size x - ultimately columns
y = cell2mat(firstHDR(6));      % image size y - ultimately rows
t=firstHDR{7};                  % Data type String
HDR = fread(fid,s(1)-50);       % get rest of header
LkUp   =     zeros(256,8,'uint8'); % Look up table for sorting out 1 bit Raw
for     i=0:255
    str     =    dec2bin(i,8);
    for     j = 1:8
         LkUp(i+1,j)=str2num(str(j));
    end
end

if strcmp(t,'U16')              % sort out data type
    A = fread(fid,[x y],'uint16');  % and read data
elseif strcmp(t,'U32')
    A = fread(fid,[x y],'uint32');
elseif  strcmp(t,'U08')
    A = fread(fid,[x y],'uint8');
elseif  strcmp(t,'R64')         % Raw, assuming 2x2 for now
    A = fread(fid,[x y],'uint16'); % Note axes 
    for k = 1:y                  % Rows
        for j = 0:(x/4-1)
            for i = 1:4
                B(j*4+i,k) = A(j*4+5-i,k); % Does 64 bit to 16 bit "endian" swap
            end
        end
    end
    if (chips == 4)         % if 4x1 then don't need to do all the stacking required for 2x2
        A = vertcat(B((x-256+1):x,1:y),B((x-2*256+1):(x-256),1:y)); % Stack chips 1 and 2, transposition coming later
        C = fliplr(flipud(B(257:512,1:y)));   % Extract chips 3 and 4, flip and invert
        D = fliplr(flipud(B(1:256,1:y)));
        D = vertcat(C,D);   % Stack chips 3 and 4
        A = horzcat(A,D);   % Now put pairs side by side
        x = 512;            % set to square image, should use layout (2x2 or Nx1) to sort out image shape
        y = 512;
    else                    % just do 1 or 4 chips
        x = 256;
        y = 256;
    end
else
    A = fread(fid,[y/8,x],'uint8'); % Read 1 bit R64 Images - this may need extra checking of bit depth at end of header
    A = A';
    ExpA = zeros(y,x);
    for i = 1:y
        for j = 1:x/8
            W = LkUp(A(i,j)+1,:);
            Col = 8*j-7;
            ExpA(i,Col:Col+7) = W;
        end
    end
    Z = fliplr(ExpA(:,1:64));
     
    for i = 65:64:x-63
        R = ExpA(:,i:i+63);
        Z = horzcat(Z,fliplr(R));
    end
    A=Z';
end
B=A';   % Transpose
% Plot if required
imagesc(B);                     % display (B)
colorbar
set(gca,'YDir','normal')        % Origin bottom left
axis equal
axis([1 x 1 y ]);
fclose(fid);
C = B;      % Return image