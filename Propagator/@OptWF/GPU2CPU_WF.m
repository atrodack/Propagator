function [ WF ] = GPU2CPU_WF( WF )
%GPU2CPU_WF Summary of this function goes here
%   Detailed explanation goes here

WF.set_field(gather(WF.field_));
WF.amp_ = gather(WF.amp_);
WF.pha_ = gather(WF.pha_);
end

