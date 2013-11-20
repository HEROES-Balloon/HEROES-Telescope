FUNCTION find_relative_source_offsets, images

reference_detector = 7

ref_image = image_remove_ring(images[*,*,1,reference_detector], center_radius = 7)
FOR idet = 0, 8-1 DO BEGIN
    IF idet NE reference_detector THEN BEGIN
        this_image = image_remove_ring(images[*,*,1,idet], center_radius = 7)
        
        res = correl_images(this_image, ref_image, xshift = 20, yshift = 20)

        tvscl, congrid(res, 300, 300)
        
        CorrMat_Analyze, res, xoffset_optimum, yoffset_optimum, max_corr
        print, xoffset_optimum, yoffset_optimum, max_corr
        ;print, size(res)
    ENDIF
ENDFOR

RETURN, -1

END