----------------------------------------------------------------
-- select consensus for each transcription factor

select
	s.site_sequence,
	bf.factor_name,
	bf.factor_ac,
	s.id,
	f.organism
from
	site s,
	site_binding_factors_expanded bf,
	factor f
where
	s.type='consensus'
and
	s.id=bf.id
and
	bf.factor_ac = f.id
and
	f.organism = 'fruit fly, Drosophila melanogaster'
;
