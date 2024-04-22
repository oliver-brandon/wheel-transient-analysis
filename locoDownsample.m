clear loco_downsample
loco_downsample = [];
for i = 1:10:107400
    x = sum(loco(i:i+10-1,:));
    loco_downsample = [loco_downsample; x];
end