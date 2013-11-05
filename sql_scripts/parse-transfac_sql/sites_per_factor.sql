----------------------------------------------------------------
-- count number of sites per factor for a selected organism
select 
	factor_ac,
	factor_name,
	count(site_binding_factors_expanded.id) site_count
from 
	site,
	site_binding_factors_expanded
where
	site.id = site_binding_factors_expanded.id
and
     site.organism = 'fruit fly, Drosophila melanogaster'
group by
	site_binding_factors_expanded.factor_ac
order by
	factor_name
--	site_count
-- into outfile 'sites_per_factor_Drosophila_melanogaster.tab'
;

