function [misalignment] = patch_transferability(source_image, target_image, source_patch, target_patch, fmaps)
    
    misalignment = sum(sum(abs((fmaps{source_image, target_image} * source_patch) - target_patch), 1));
    misalignment = misalignment + sum(sum(abs((fmaps{target_image, source_image} * target_patch) - source_patch), 1));

end
