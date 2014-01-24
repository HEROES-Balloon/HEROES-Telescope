FUNCTION convert_detcoord_to_physicalcoord, rawx, rawy, detector_number

;PURPOSE: Converts from detectors (raw) coordinates to physical coordinates.
;
;WRITTEN: Steven Christe (14-Jan-2014)

detector_raw_to_arcsec = 3.4

detector_raw_center = fltarr(2, 8)
detector_raw_center[*, 0] = [320., 280.]
detector_raw_center[*, 1] = [320., 300.]
detector_raw_center[*, 2] = [295., 250.]
detector_raw_center[*, 3] = [330., 300.]
detector_raw_center[*, 4] = [280., 280.]
detector_raw_center[*, 5] = [305., 280.]
detector_raw_center[*, 6] = [280., 290.]
detector_raw_center[*, 7] = [330., 290.]

x_arcsec = (rawx - detector_raw_center[0]+0.5) * detector_raw_to_arcsec
y_arcsec = (rawy - detector_raw_center[1]+0.5) * detector_raw_to_arcsec

RETURN, [x_arcsec, y_arcsec]

END