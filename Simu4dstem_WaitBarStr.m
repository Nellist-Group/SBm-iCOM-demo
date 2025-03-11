classdef Simu4dstem_WaitBarStr
    %WAITBARSTR 此处显示有关此类的摘�?
    %   version: 0.1.1
    
    properties
        NowValue = 0;
        DoneUnit = '#';
        BlankUnit = '~';
        MaxLength = 50;
        Commonts = [];
    end
    
    methods
        function Init(self,Comm)
            self.NowValue = 0;
            self.Commonts = Comm;
            fprintf('%s\t@ %s\n',Comm,datestr(now));
            fprintf('  0%%');
            fprintf('|%s|',self.DuplicateUnit(self.BlankUnit,self.MaxLength));
        end
        
        function Update(self,R)
            % R between 0 and 1
            self.NowValue = R;
            if R < 1
                UnitNum = floor(R*self.MaxLength);
            elseif R >= 1
                UnitNum = self.MaxLength;
            end
            self.DeleteLine(self.MaxLength+6);
            NowStr = [self.DuplicateUnit(self.DoneUnit,UnitNum),...
                self.DuplicateUnit(self.BlankUnit,self.MaxLength-UnitNum)];
            fprintf('% 3u%%',floor(self.NowValue*100));
            fprintf('|%s|',NowStr);
            if R >= 1
                fprintf('\nFinished @ %s\n',datestr(now));
            end
        end
        
        function DUnit = DuplicateUnit(self,Unit,Num)
            DUnit = [];
            for ii = 1:Num
                DUnit = strcat(DUnit,Unit);
            end
        end
        
        function DeleteLine(self,Num)
            for ii = 1:Num
                fprintf('\b');
            end
        end
    end
    
end 

