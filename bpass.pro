PRO BPass,hora,smoooth,result

delta = 6.          ; quantidade de pontos num intervalo de uma hora
pcl = 3*delta    	; periodos abaixo de pcl passarao (em horas)
pch = (3/6.)*delta	; periodos acima de pch passarao

pl = float(pcl*delta/delta)		;pl and ph are the appropriate cutoff periods in time units
ph = float(pch*delta/delta)

in_fl = (1/pl)	 ;Cut-in frequency (i.e. filter has zero response from
;                 f=[0,fc1]).  In units of cycles per data interval.
;                 Floating/double scalar.  Not changed by procedure.

in_fh = (1/ph)	 ;Cut-out frequency (i.e. filter has zero response from
;                 f=[fc2,Nyquist]).  In units of cycles per data interval.
;                 Floating/double scalar.  Not changed by procedure.

n = FIX(1.3/(in_fh - in_fl))

in_n = 4 * n         ;The parameter n, which defines how many 2n+1 weights to
;                  use in the filter (so weights go from k=-n to n).

;print, in_fl, ' - ', in_fh

Result = LANCZOS_BANDPASS(smoooth, in_fl, in_fh, in_n)

END