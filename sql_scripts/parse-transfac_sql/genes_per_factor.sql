----------------------------------------------------------------
-- target genes for each factor
select distinct 
     site.gene_name,
     site_binding_factors_expanded.factor_name,
     site.gene_ac,
     site_binding_factors_expanded.factor_ac
from 
     site,
     site_binding_factors_expanded
where 
     site.organism = 'fruit fly, Drosophila melanogaster'
and
	site.id = site_binding_factors_expanded.id
order by
	site_binding_factors_expanded.factor_name
into outfile 'gene_factor_Drosophila_melanogaster.tab'
;

----------------------------------------------------------------
-- Count target genes for each factor
select 
     site_binding_factors_expanded.factor_name,
     site_binding_factors_expanded.factor_ac,
     count(distinct site.gene_ac) gene_count
from 
     site,
     site_binding_factors_expanded
where 
     site.organism = 'fruit fly, Drosophila melanogaster'
and
	site.id = site_binding_factors_expanded.id
group by 
	factor_ac
order by
	gene_count
--into outfile 'genes_per_factor_Drosophila_melanogaster.tab'
;
