;+
; NAME:
;   dcomdisdz
; PURPOSE:
;   Compute differential comoving line-of-sight distances (for c/H_0=1).
; CALLING SEQUENCE:
;   dDdz= dcomdisdz(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   differential comoving distance DD/dz in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;-
function dcomdisdz, z,OmegaM,OmegaL
  return, (1.0/sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL))
end
