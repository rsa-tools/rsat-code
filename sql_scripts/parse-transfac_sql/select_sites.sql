----------------------------------------------------------------
-- select sites for a selected organism
-- export columns in the appropriate order for feature-map
select 
     site.gene_name,
     'site',
     site_binding_factors_expanded.factor_name,
     'DR',
     site.position_first,
     site.position_last,
     site.site_sequence,
     '',
     site.id,
     site.gene_ac,
     site_binding_factors_expanded.factor_ac,
     site.organism
from 
     site,
     site_binding_factors_expanded
where 
     site.organism = 'fruit fly, Drosophila melanogaster'
and
	site.id = site_binding_factors_expanded.id
into outfile
      'transfac_sites_Drosophila_melanogaster.tab'
--      '/Users/jvanheld/rsa-tools/public_html/tmp/transfac_sites_Drosophila_melanogaster.tab'
;


