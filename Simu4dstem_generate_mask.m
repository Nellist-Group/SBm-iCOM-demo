function Mask = Simu4dstem_generate_mask(MaskShape,PxlNum,varargin)
% generateMask1 --> center in ceil((size+1)/2)
% MaskShape:'center-dot'/
%   if MaskShape = 'center-dot'
%       varargin{1} = r in pxl
%   elseif MaskShape = 'center-ring'
%       varargin{1} = inner r inpxl
%       varargin{2} = outer r inpxl
%   elseif MaskShape = 'center-dark-dot'
%       varargin{1} = r in pxl
%   elseif MaskShape = 'ramp'
%       no varargin
%   elseif MaskShape = 'center-ring-sector'
%       varargin{1} = inner r in pxl
%       varargin{2} = outer r in pxl
%       varargin{3} = start angle in [-pi, pi]
%       varargin{4} = angle range in rad, positive or negtive

tx = -(floor(PxlNum(1)/2)):(floor((PxlNum(1)-1)/2));
ty = -(floor(PxlNum(2)/2)):(floor((PxlNum(2)-1)/2));
[y,x] = meshgrid(ty,tx);
switch MaskShape
    case 'center-dot'
        Mask = zeros(PxlNum(1),PxlNum(2));
        Mask((x.^2+y.^2) <= varargin{1}^2) = 1;
    case 'center-ring'
        Mask = ones(PxlNum(1),PxlNum(2));
        Mask((x.^2+y.^2) <= varargin{1}^2) = 0;
        Mask((x.^2+y.^2) > varargin{2}^2) = 0;
    case 'center-dark-dot'
        Mask = ones(PxlNum(1),PxlNum(2));
        Mask((x.^2+y.^2) <= varargin{1}^2) = 0;
    case 'ramp'
        top = sqrt(PxlNum(1).^2+PxlNum(2).^2) / 2;
        Mask = -sqrt(x.^2+y.^2)/top + 1;
    case 'center-ring-sector'
        Mask = ones(PxlNum(1),PxlNum(2));
        Mask((x.^2+y.^2) <= varargin{1}^2) = 0;
        Mask((x.^2+y.^2) > varargin{2}^2) = 0;
        ti = y + 1j .* x;
        ta = angle(ti);
        EndAngle = varargin{3} + varargin{4};
        if EndAngle <= pi && EndAngle >= -pi
            if varargin{4} > 0
                AngleMask = (ta > varargin{3} & ta < EndAngle);
            else
                AngleMask = (ta < varargin{3} & ta > EndAngle);
            end
        elseif EndAngle > pi
            AngleMask = (ta > varargin{3} | ta + 2*pi < EndAngle);
        elseif EndAngle < -pi
            AngleMask = (ta < varargin{3} | ta - 2*pi > EndAngle);
        end
        Mask = Mask .* AngleMask;
    otherwise
        error('wrong shape input');
end




end

