----------------------------------------------------------------
-- select sites for a selected organism
-- export columns in the appropriate order for feature-map
select 
     s.gene_name,
     'site',
     bf.factor_name,
     'DR',
     s.position_first,
     s.position_last,
     s.site_sequence,
     '',
     s.id,
     s.gene_ac,
     bf.factor_ac,
     s.organism
from 
     site s,
     site_binding_factors_expanded bf
where 
     s.organism = 'fruit fly, Drosophila melanogaster'
and  s.id = bf.id
-- and  bf.factor_name like '%Prd%'
into outfile
     'transfac_sites_Drosophila_melanogaster.tab'
--     '/Users/jvanheld/rsa-tools/public_html/tmp/transfac_sites_Drosophila_melanogaster.tab'
;


