function record = make_motion_trial(time,handXYZ,jointAngles,jointNames,sync)
%MAKE_MOTION_TRIAL Explicit example adapter: hand in mm, joint angles in rad.
% Change the conversions/frame names to match verified acquisition metadata.
% record is assigned to MotionData.P4.trial2, for example.
if nargin<5, sync=struct(); end
if size(handXYZ,1)~=numel(time)||size(handXYZ,2)~=3
    error('example:Shape','handXYZ must be N-by-3 signed coordinates.');
end
blank=struct('name','','kind','','unit','','values',[],'frame','','valid',[], ...
    'noiseSD',[],'limitRMS',NaN,'calibrationVerified',false,'calibrationPolicyHash','', ...
    'reviewed',[],'muscles',1:15); % Empty target flag inherits trial.reviewed.
targets=repmat(blank,1,1+numel(jointNames));
targets(1).name='hand_torso'; targets(1).kind='translation'; targets(1).unit='mm';
targets(1).values=handXYZ; targets(1).frame='thorax'; targets(1).valid=all(isfinite(handXYZ),2);
if size(jointAngles,1)~=numel(time)||size(jointAngles,2)~=numel(jointNames)
    error('example:Shape','jointAngles is N-by-number_of_jointNames.');
end
for j=1:numel(jointNames)
    targets(j+1).name=jointNames{j}; targets(j+1).kind='joint_angle'; targets(j+1).unit='rad';
    targets(j+1).values=jointAngles(:,j); targets(j+1).frame='calibrated_anatomical_joint';
    targets(j+1).valid=isfinite(jointAngles(:,j));
end
record=struct('time',time(:),'targets',targets,'sync',sync,'reviewed',false);
end
