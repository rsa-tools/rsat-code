----------------------------------------------------------------
-- Count number of sites per organism
-- 
select 
     organism,
     count(id) sites
from 
     site
group by 
      organism
order by
      sites
;
