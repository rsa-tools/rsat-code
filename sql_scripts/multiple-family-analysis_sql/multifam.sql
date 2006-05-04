----------------------------------------------------------------
--
-- Process results of multiple-family-analysis
--
----------------------------------------------------------------

----------------------------------------------------------------
-- number of families per condition
select distinct 
       suffix,
       count(id)
from
	family
group by
      suffix
;

----------------------------------------------------------------
-- Count genes per family
--
select 
       family.suffix suffix,
       family.id family,
       count(family_genes.genes) genenb
from 
     family,
     family_genes
where
	family.id = family_genes.id
group by 
      suffix,
      family
order by
      suffix,
      genenb
;

