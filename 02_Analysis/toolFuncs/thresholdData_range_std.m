function [bad_segments,dead_segments_perc,noise_segments_perc,dead_segments,noise_segments] = thresholdData_range_std(data,srate,cfg_thresh)

% data must be organizes channelsxsamples
if nargin<3
    cfg_thresh.segments_size_sec   = 2;
    cfg_thresh.dead_thresh_var     = 1;
    cfg_thresh.variance_thresh     = 250;
    cfg_thresh.range_thresh        = 500;
end

Y = movstd(data,srate*cfg_thresh.segments_size_sec,0,2,'Endpoints',"fill");

%
seg_compl           = ones(1,srate*cfg_thresh.segments_size_sec);
% %deadsegments
dead_segments       = Y<cfg_thresh.dead_thresh_var;
lr                  = conv2(dead_segments,seg_compl,'full');
rl                  = conv2(fliplr(dead_segments),seg_compl,'full');
dead_segments       = (lr(:,1:size(dead_segments,2))+fliplr(rl(:,1:size(dead_segments,2))))>0;
dead_segments_perc  = sum(dead_segments,2)/size(Y,2);

% horrible noise semgemtns
noise_segments      = Y>cfg_thresh.variance_thresh | abs(data)>cfg_thresh.range_thresh;
lr                  = conv2(noise_segments,seg_compl,'full');
rl                  = conv2(fliplr(noise_segments),seg_compl,'full');
noise_segments       = (lr(:,1:size(noise_segments,2))+fliplr(rl(:,1:size(noise_segments,2))))>0;

noise_segments_perc = sum(noise_segments,2)./size(Y,2);

bad_segments = any(dead_segments+noise_segments);

starts = find(diff(bad_segments)==1);       % this gives one sample before bad segment
ends   = find(diff(bad_segments)==-1)+1;    % this gives one sample after bad segment

if bad_segments(1)==1
    starts = [1,starts];
end

if bad_segments(end)==1
    ends = [ends,length(bad_segments)];
end

bad_segments = [starts',ends'];
