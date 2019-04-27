;+
; NAME:
;   angdidis
; PURPOSE:
;   Compute angular diameter distancea (for c/H_0=1).
; CALLING SEQUENCE:
;   D= angdidis(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   angular diameter distance in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   propmotdis()
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;-
function angdidis, z,OmegaM,OmegaL
  return, propmotdis(z,OmegaM,OmegaL)/(1.0+z)
end
