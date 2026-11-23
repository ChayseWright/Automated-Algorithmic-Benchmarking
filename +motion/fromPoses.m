function record = fromPoses(time,segments,referenceName,handName,unit,jointPairs)
%FROMPOSES Explicit 6-DOF adapter. Quaternions are [w x y z], active rotations.
% Relative rotation vectors are NOT calibrated anatomical joint DOFs.
if nargin<6, jointPairs=cell(0,2); end
ref=find(strcmp({segments.name},referenceName),1); hand=find(strcmp({segments.name},handName),1);
if isempty(ref)||isempty(hand), error('tbench:MotionPose','Reference/hand segment name is missing.'); end
N=numel(time); poses=segments;
for k=1:numel(poses)
    if ~isequal(size(poses(k).position),[N 3])||~isequal(size(poses(k).quaternion),[N 4])
        error('tbench:MotionPose','Positions are N-by-3 and quaternions N-by-4.');
    end
    q=double(poses(k).quaternion); normq=sqrt(sum(q.^2,2)); good=all(isfinite(q),2)&normq>1e-12;
    good=good&all(isfinite(poses(k).position),2);
    if isfield(poses(k),'valid')&&~isempty(poses(k).valid)
        mask=poses(k).valid;
        if numel(mask)~=N||any(~isfinite(mask(:)))||any(mask(:)~=0&mask(:)~=1)
            error('tbench:MotionPose','Segment valid must contain N finite 0/1 flags.');
        end
        good=good&logical(mask(:));
    end
    q=bsxfun(@rdivide,q,normq); q(~good,:)=NaN; poses(k).quaternion=q;
    poses(k).position(~good,:)=NaN;
end
q=conjugate(poses(ref).quaternion);
translation=rotate(q,poses(hand).position-poses(ref).position);
targets=makeTarget([handName '_relative_' referenceName],'translation',unit,translation,referenceName);
for k=1:size(jointPairs,1)
    parent=find(strcmp({poses.name},jointPairs{k,1}),1); child=find(strcmp({poses.name},jointPairs{k,2}),1);
    if isempty(parent)||isempty(child), error('tbench:MotionPose','A joint-pair segment is missing.'); end
    relative=product(conjugate(poses(parent).quaternion),poses(child).quaternion);
    negative=relative(:,1)<0; relative(negative,:)=-relative(negative,:);
    v=relative(:,2:4); s=sqrt(sum(v.^2,2)); angle=2*atan2(s,relative(:,1)); factor=2*ones(size(s));
    good=s>1e-12; factor(good)=angle(good)./s(good); values=bsxfun(@times,v,factor);
    values(angle>pi-.02,:)=NaN; % The logarithm chart is singular near pi.
    targets(end+1)=makeTarget([jointPairs{k,2} '_relative_' jointPairs{k,1} '_rotation_vector'], ...
        'rotation','rad',values,jointPairs{k,1}); %#ok<AGROW>
end
record=struct('time',time(:),'targets',targets,'sync',struct(),'reviewed',false);
end
function t=makeTarget(name,kind,unit,v,frame)
t=struct('name',name,'kind',kind,'unit',unit,'values',v,'frame',frame,'valid',all(isfinite(v),2));
end
function q=conjugate(q)
q(:,2:4)=-q(:,2:4);
end
function c=product(a,b)
c=[a(:,1).*b(:,1)-sum(a(:,2:4).*b(:,2:4),2), ...
    bsxfun(@times,a(:,1),b(:,2:4))+bsxfun(@times,b(:,1),a(:,2:4))+cross(a(:,2:4),b(:,2:4),2)];
end
function v=rotate(q,v)
u=q(:,2:4); v=v+2*bsxfun(@times,q(:,1),cross(u,v,2))+2*cross(u,cross(u,v,2),2);
end
