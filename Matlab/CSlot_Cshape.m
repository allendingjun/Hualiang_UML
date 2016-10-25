
function CSlot_Cshape(fid,r1,alp1,theta1,w1,r,alp,theta,w,rc,thickness,P,NUM)

%%% r1,alp1,theta1,w1: for outer C-Slot
%%% r,alp,theta,w: for inner C-shape
%%% NUM is the name of the created cell

% Cshape
CshapeName=['Cshape' num2str(NUM)];
hfssRectangle(fid, CshapeName, 'X', [0,r-w,0], w, thickness, 'um'); %%% creat a rectangle: thickness is Z direction length
hfssSweepAroundAxis(fid, CshapeName, 'Z', 360-alp);
hfssRotate(fid, {CshapeName,'Rec3'}, 'Z', theta-alp/2);

% CSlot
CSlotName=['Cslot' num2str(NUM)];
hfssRectangle(fid, CSlotName, 'X', [0,r1-w1,0], w1, thickness, 'um'); %%% creat a rectangle: thickness is Z direction length
hfssSweepAroundAxis(fid, CSlotName, 'Z', 360-alp1);
hfssRotate(fid, {CSlotName,'Rec4'}, 'Z', theta1-alp1/2);

% circle to be cut
CutCircleName=['CutCircle' num2str(NUM)];
hfssCylinder(fid, CutCircleName, 'Z', [0,0,0], rc, thickness, 'um');

%% The ground
GroundName=['Ground' num2str(NUM)];
hfssBox(fid, GroundName, [-P/2, -P/2, 0],[P, P, thickness], 'um');

hfssSubtract(fid,GroundName,CutCircleName);
hfssSubtract(fid,GroundName,CSlotName);

hfssUnite(fid, {GroundName, CshapeName});



