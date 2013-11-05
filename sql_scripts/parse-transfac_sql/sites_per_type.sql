----------------------------------------------------------------
-- Count number of sites per organism
-- 
select 
     type,
     count(id) sites
from 
     site
group by 
      type
order by
      sites
;

----------------------------------------------------------------
-- Select sites with no organism, and not artificial
-- 
select 
     type,
     count(id) sites
from 
     site
where
	organism = '<NULL>'
group by 
      type
order by
      sites
;


